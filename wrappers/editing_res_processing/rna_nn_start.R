# 
# if(!require(LncFinder)) {install.packages("LncFinder", repos = getCRANmirrors()$URL[1])}
# library("LncFinder")

library(data.table,verbose = F)
library(rtracklayer,verbose = F)
library(stringr,verbose = F)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(seqinr)
library(tictoc,verbose = F)
library(parallel,verbose = F)






read_single_chr_bam_as_sam <- function(chr_name,bam_file){

   sam <- fread(cmd = paste0("samtools view ", bam_file," ",chr_name," | awk -v OFS='\t' '$0 ~ /MD:Z:[0-9]*[A-Z]/{match($0,/MD:Z:[^\t]+/); print $1, $2, $3, $4, $6, $10, $11, substr($0,RSTART,RLENGTH)}'"))
   
   # #to_test
   # sam <- fread("/mnt/ssd/ssd_1/workspace/vojta/A2I_edit_EXO10_KD/preprocessed_u87_adar_kd_rep2.tsv")

   
   setnames(sam,c("qname", "flag", "chr", "start", "cigar", "read", "mapq", "MD"))
   sam[,end := start + nchar(read)]
   sam[,chr := as.character(chr)]
   
   return(sam)
}


gene_annotation_and_multimapped_split <- function(sam,annot_ref){

   sam[,row_id := seq_along(chr)]
   sam <- foverlaps(sam,annot_ref)
   setnames(sam,c("start","end","i.start","i.end"),c("gene_start","gene_end","mapping_start","mapping_end"))
   sam[is.na(gene_id),gene_id := "ncRNA"]
   sam[gene_id == "ncRNA",strand := "+"]
   sam[,per_gene_read_count := .N,by = c("gene_id")]
   setorder(sam,-per_gene_read_count) 
   sam <- unique(sam,by = c("row_id","strand"))
   sam[,per_gene_read_count := NULL]
   sam[is.na(gene_id),gene_id := "ncRNA"]
   
   #separate genes by strand
   strand_cast_tab <- dcast.data.table(sam,formula = row_id ~ strand,value.var = "gene_id",fill = "ncRNA")
   setnames(strand_cast_tab,c("row_id","fwd_gene_id","rev_gene_id"))
   sam[,gene_id := NULL]
   sam[,strand := NULL]
   sam <- unique(sam,by = c("row_id"))
   setkey(sam,row_id)
   sam <- strand_cast_tab[sam]

   return(sam)

}





create_cigar_table <- function(sam){
   #getting cigar information
   # Op BAM Description Consumes
   # query
   # Consumes
   # reference
   # M 0 alignment match (can be a sequence match or mismatch) yes yes
   # I 1 insertion to the reference yes no
   # D 2 deletion from the reference no yes
   # N 3 skipped region from the reference no yes
   # S 4 soft clipping (clipped sequences present in SEQ) yes no
   # H 5 hard clipping (clipped sequences NOT present in SEQ) no no
   # P 6 padding (silent deletion from padded reference) no no
   # = 7 sequence match yes yes
   # X 8 sequence mismatch yes yes
   
   
   cigar_tab <- str_match_all(sam$cigar, "([0-9]+)([SMIND]+)")
   cigar_tab <- data.table(row_id = rep(sam$row_id,sapply(cigar_tab, nrow)),
                           length = as.numeric(unlist(lapply(cigar_tab,function(x) x[,2]))),
                           type = unlist(lapply(cigar_tab,function(x) x[,3])))
   cigar_tab[,cons_que := F]
   cigar_tab[,cons_ref := F]
   cigar_tab[type %in% c("M","I","S"),cons_que := T]
   cigar_tab[type %in% c("M","D","N"),cons_ref := T]
   cigar_tab[,cons_que_length := length]
   cigar_tab[cons_que == F,cons_que_length := 0]
   cigar_tab[,que_start_pos := cumsum(cons_que_length) - cons_que_length +1,by = row_id]
   setkey(cigar_tab,row_id)
   
   
   return(cigar_tab)
}
   
   

create_mismatch_table <- function(cigar_tab,sam){
      
   #MD_information
   sam[, MD := gsub("MD:Z:", "", MD)]
   sam[grepl("^[ACTG]",MD),MD := paste0('0', MD)]
   str_match <- str_match_all(sam$MD, "([0-9]+)([A-Z]|\\^[A-Z]+)")
   
   mismatch_tab <- data.table(row_id = rep(sam$row_id,sapply(str_match, nrow)),
                               length = as.numeric(unlist(lapply(str_match,function(x) x[,2]))),
                               ref = unlist(lapply(str_match,function(x) x[,3])))
   
   mismatch_tab[,nchar_ref := nchar(ref)]
   mismatch_tab[nchar_ref > 1,char_size_offset := 0]
   mismatch_tab[nchar_ref == 1,char_size_offset := 1]
   mismatch_tab[,nchar_ref := NULL]
   mismatch_tab[,read_pos := cumsum(length + char_size_offset),by = row_id]
   mismatch_tab[,char_size_offset := NULL]
   mismatch_tab[,ref_pos := read_pos + sam[J(mismatch_tab$row_id)]$mapping_start - 1]
   
   mismatch_tab <- mismatch_tab[nchar(ref) == 1]
   
   setkey(mismatch_tab,row_id,read_pos)
   mismatch_tab <- merge(mismatch_tab,cigar_tab[cons_ref == F,.(row_id,consume_length = length,cigar_pos = que_start_pos)],by = "row_id",all.x = T,allow.cartesian=TRUE)
   
   mismatch_tab[is.na(consume_length),consume_length := 0]
   mismatch_tab[is.na(cigar_pos),cigar_pos := 0]
   mismatch_tab <- mismatch_tab[,.(ref = ref[1],ref_pos = ref_pos[1],read_pos_new = read_pos[1] + sum(consume_length[cigar_pos < read_pos])),by = .(row_id,read_pos)]
   mismatch_tab[,read_pos := read_pos_new]
   mismatch_tab[,read_pos_new := NULL]
   
   mismatch_tab[,alt := stringr::str_sub( sam[J(mismatch_tab$row_id)]$read,read_pos,read_pos)]
   mismatch_tab[,base_quality := stringr::str_sub( sam[J(mismatch_tab$row_id)]$mapq,read_pos,read_pos)]
   
   mismatch_tab <- merge(mismatch_tab,cigar_tab[cons_que == F,.(row_id,consume_length = length,cigar_pos = que_start_pos)],by = "row_id",all.x = T,allow.cartesian=TRUE)
   mismatch_tab[is.na(consume_length),consume_length := 0]
   mismatch_tab[is.na(cigar_pos),cigar_pos := 0]
   
   mismatch_tab <- mismatch_tab[,.(read_pos = read_pos[1],ref = ref[1],ref_pos_new = ref_pos[1] + sum(consume_length[cigar_pos <= read_pos])),by = .(row_id,ref_pos)]
   mismatch_tab[,ref_pos := ref_pos_new]
   mismatch_tab[,ref_pos_new := NULL]
   
   mismatch_tab[,alt := stringr::str_sub( sam[J(mismatch_tab$row_id)]$read,read_pos,read_pos)]
   mismatch_tab[,base_quality := stringr::str_sub( sam[J(mismatch_tab$row_id)]$mapq,read_pos,read_pos)]
   mismatch_tab <- mismatch_tab[ref != "N",]
   mismatch_tab <- mismatch_tab[alt != "N",]
   mismatch_tab <- mismatch_tab[alt != ref,]
   
   mismatch_tab <- sam[,.(row_id,chr,fwd_gene_id,rev_gene_id,mapping_start,mapping_end)][mismatch_tab]
   mismatch_tab[,per_map_mismatch_count := .N,by = row_id]
   
   
   return(mismatch_tab)
   
   #test code
   # test_tab <- dt_AG
   # tab_ranges <- GRanges(paste0("chr",test_tab$chr),
   #                       IRanges(start=test_tab$ref_pos, end = test_tab$ref_pos),
   #                       strand="+")
   # 
   # tab_seq <- as.data.table(Views(BSgenome.Hsapiens.UCSC.hg38,tab_ranges))
   # 
   # mismatch_tab <- merge(mismatch_tab,tab_seq[,.(chr = gsub("chr","",seqnames),ref_pos = start,orig_ref = dna)],by = c("chr","ref_pos"))
   # bad_tab <- mismatch_tab[ref != orig_ref]
   # bad_tab[,n := .N,by = pos]
   # bad_tab[n == max(bad_tab$n)]
   
   
}

filter_editing_specific_variants <- function(mismatch_tab,edit_setings = "+A>G;-T>C|+G>A;-C>T"){
   processed_edit_setings <- strsplit(edit_setings,split = "\\|")[[1]]
   processed_edit_setings <- unlist(strsplit(processed_edit_setings,split = ";"))
   all_filter_texts <- sapply(processed_edit_setings,function(text){
      text_vec <- strsplit(text,split = "")[[1]]
      if(text_vec[1] == "+"){
         return(paste0("mismatch_tab[ref == '",text_vec[2],"' & alt == '",text_vec[4],"',gene_id := fwd_gene_id]"))
      } else {
         return(paste0("mismatch_tab[ref == '",text_vec[2],"' & alt == '",text_vec[4],"',gene_id := rev_gene_id]"))
      }
   })
   for(filter_text in all_filter_texts){
      eval(parse(text = filter_text))
   }
   mismatch_tab <- mismatch_tab[!is.na(gene_id)]
   mismatch_tab[,fwd_gene_id := NULL]
   mismatch_tab[,rev_gene_id := NULL]
   
   return(mismatch_tab)
}

base_qualiy_score <- function(mismatch_tab){
   # adding 'probability score' of real editing):   
   mismatch_tab[, Q := as.numeric(sapply(base_quality, function(x) as.integer(charToRaw(as.character(x)))-33))]
   mismatch_tab[, P_error := 10^-(Q/10)]
   mismatch_tab[, qual_prob := (1-P_error) * 100]
   mismatch_tab[,c("base_quality","Q","P_error") := NULL]
   
   return(mismatch_tab)
}

known_editing_sites <- function(mismatch_tab,chrom_known_edit_tab){
   #Finding known editing sites
   chrom_known_edit_tab[, known_editting := TRUE]
   mismatch_tab <- chrom_known_edit_tab[mismatch_tab]
   mismatch_tab[is.na(known_editting),known_editting := F]

   return(mismatch_tab)
}


known_snps <- function(mismatch_tab,chrom_snp_tab){
   ## check SNPs
   mismatch_tab <- chrom_snp_tab[mismatch_tab]
   mismatch_tab[is.na(pop_AF),pop_AF := 0]
   
   return(mismatch_tab)
}

single_var_tab_cluster_computation <- function(cluster_tab,strand,var_dist_window){
   if(strand == "+"){
      cluster_tab[,gene_id := fwd_gene_id]
   } else {
      cluster_tab[,gene_id := rev_gene_id]
   }
   cluster_tab[,diff_vec := c(diff(sort(cluster_tab$ref_pos)),Inf)]
   diff_vec <- diff(sort(cluster_tab$ref_pos))
   diff_vec[diff_vec < var_dist_window] <- 0
   diff_vec[diff_vec > 0] <- seq_along(diff_vec[diff_vec > 0])
   rle_diff_vec <- rle(diff_vec)
   cluster_tab[,cluster := c(1,rep(seq_along(rle_diff_vec$lengths),times = rle_diff_vec$lengths))]
   cluster_tab[cluster %in% which(rle_diff_vec$values > 0),cluster := cluster + 1]
   cluster_tab[,gene_cluster_size := .N ,by = .(gene_id,cluster)]
   cluster_tab[,is_correct_gene := gene_id == gene_id[which.max(gene_cluster_size)],by = cluster]
   cluster_tab[,gene_cluster_size := NULL]
   cluster_tab <- cluster_tab[is_correct_gene == T]
   cluster_tab[,is_correct_gene := NULL]
   cluster_tab[,var_count := .N,by = cluster]
   
   
   cluster_tab <- cluster_tab[var_count > 1]
   cluster_tab <- cluster_tab[,.(chr = chr[1],
                                 var_count = var_count[1],
                                 cluster_size = max(ref_pos) - min(ref_pos),
                                 cluster_start_pos = min(ref_pos),
                                 cluster_end_pos = max(ref_pos),
                                 mean_median_bias = mean(table(ref_pos)) - median(table(ref_pos)),
                                 position_count = length(unique(ref_pos)),
                                 gene_id = gene_id[1]),by = cluster]
   cluster_tab[,cluster := NULL]
   return(cluster_tab)
}


process_single_chromosome <- function(bam_file,chrom_ref,chrom_name,chrom_known_edit_tab,chrom_snp_tab,edit_variation,var_dist_window,cluster_window_extension){
   
   sam <- read_single_chr_bam_as_sam(chrom_name,bam_file) 
   sam <- gene_annotation_and_multimapped_split(sam,chrom_ref[,.(chr,start,end,gene_id,strand)])
   
   cigar_tab <- create_cigar_table(sam)
   
   # coverage_tab <- cigar_tab[cons_ref == T][sam[,.(row_id,chr,mapping_start)]]
   # coverage_tab[,mapping_start := mapping_start + c(0,cumsum(head(length,-1))),by = row_id]
   # coverage_tab <- coverage_tab[type == "M",.(chr,mapping_start,mapping_end = mapping_start + length)]
   # setkey(coverage_tab)
   
   mismatch_tab <- create_mismatch_table(cigar_tab,sam)
   
   setkey(mismatch_tab,chr,ref_pos)
   mismatch_tab <- base_qualiy_score(mismatch_tab)
   mismatch_tab <- known_editing_sites(mismatch_tab,chrom_known_edit_tab)
   mismatch_tab <- known_snps(mismatch_tab,chrom_snp_tab)
   
   #quality and population filter
   filter_mismatch_tab <- mismatch_tab[qual_prob > 96]    
   filter_mismatch_tab <- filter_mismatch_tab[pop_AF == 0]
   
   edit_signal_vars <- strsplit(edit_variation,split = ";")[[1]]
   edit_signal_vars <- matrix(unlist(tstrsplit(edit_signal_vars,split = "")[c(1,2,4)]),nrow = 2)
   
   edit_cluster <- apply(edit_signal_vars,1,function(var_line) single_var_tab_cluster_computation(filter_mismatch_tab[ref == var_line[2] & alt == var_line[3]],strand = var_line[1],var_dist_window) )
   names(edit_cluster) <- apply(edit_signal_vars,1,paste,collapse = "")
   edit_cluster <- rbindlist(edit_cluster,use.names = T,idcol = "variant")
   
   edit_cluster[,extend_cluster_start_pos := cluster_start_pos - cluster_window_extension]
   edit_cluster[,extend_cluster_end_pos := cluster_end_pos + cluster_window_extension]
   
   setkey(edit_cluster,chr,extend_cluster_start_pos,extend_cluster_end_pos)
   
   # mismatch_cluster_overlap <- foverlaps(edit_cluster,setkey(mismatch_tab[,.(chr,cluster_start_pos = ref_pos,cluster_end_pos = ref_pos,ref,alt)],chr,cluster_start_pos,cluster_end_pos))
   mismatch_cluster_overlap <- foverlaps(edit_cluster,
                                         setkey(mismatch_tab[,.(chr,extend_cluster_start_pos = ref_pos,extend_cluster_end_pos = ref_pos,ref,alt)],chr,extend_cluster_start_pos,extend_cluster_end_pos))
   mismatch_cluster_overlap[,extend_cluster_end_pos := NULL]
   setnames(mismatch_cluster_overlap,c("extend_cluster_start_pos","i.extend_cluster_start_pos","i.extend_cluster_end_pos"),c("ref_pos","extend_cluster_start_pos","extend_cluster_end_pos"))
   mismatch_cluster_overlap[,not_extended_pos := ref_pos >= extend_cluster_start_pos + cluster_window_extension & ref_pos <= extend_cluster_end_pos - cluster_window_extension]
   mismatch_cluster_overlap <- mismatch_cluster_overlap[,.(all_vars = sum(not_extended_pos),
                                                           all_vars_pos = length(unique(ref_pos[not_extended_pos])),
                                                           all_vars_extended = .N,
                                                           all_vars_pos_extended = length(unique(ref_pos))),by = .(chr,extend_cluster_start_pos,extend_cluster_end_pos)]
   
   
   edit_cluster <- edit_cluster[mismatch_cluster_overlap]
   
   # setkey(edit_cluster,chr,cluster_start_pos,cluster_end_pos)
   # 
   # mismatch_coverage_overlap <- foverlaps(edit_cluster,coverage_tab,by.x=c("chr","cluster_start_pos","cluster_end_pos"), by.y=c("chr","mapping_start","mapping_end"))
   # 
   # mismatch_coverage_overlap[,overlap := cluster_size - pmax.int(0,mapping_start - cluster_start_pos) - pmax.int(0,cluster_end_pos - mapping_end)]
   # edit_cluster <- edit_cluster[mismatch_coverage_overlap[,.(read_basepairs = sum(overlap)),by = .(chr,cluster_start_pos,cluster_end_pos)]]
   # 
   # gene_coverage_tab <- foverlaps(coverage_tab,chrom_ref[,.(chr,start,end,gene_id)],by.x=c("chr","mapping_start","mapping_end"), by.y=c("chr","start","end"))
   # gene_coverage_tab[,overlap := mapping_end - mapping_start]
   # gene_coverage_tab[!is.na(gene_id), overlap := overlap - pmax.int(0,start - mapping_start) - pmax.int(0,mapping_end - end)]
   # gene_coverage_tab[is.na(gene_id), gene_id := "ncRNA"]
   # gene_coverage_tab <- gene_coverage_tab[,.(size = end[1] - start[1],read_basepairs = sum(overlap)),by = gene_id]
   
   
   return(list(edit_cluster = edit_cluster,mismatch_tab = filter_mismatch_tab))
   
}

run_all <- function(args){
   bam_file <- args[1]
   gtf_file <- args[2]
   known_edit_file <- args[3]
   snp_file <- args[4]
   output_file_prefix <- args[5]
   # edit_setings = "+A>G;-T>C|+G>A;-C>T"
   edit_setings = "+A>G;-T>C|+T>A;-C>G"
   edit_variation = "+A>G;-T>C"
   var_dist_window <- 50
   
   min_gene_vars = 20
   smooth_coef = 2
   cluster_window_extension = 20
   
   tic("full_sample_analysis")
   
   #read gtf_file
   feat_type <- "gene"
   annotate_by<- c("gene_name","gene_id", "strand")
   ref <- as.data.table(rtracklayer::import(gtf_file, feature.type = feat_type))[, c("seqnames","start","end",annotate_by), with=F]
   setnames(ref,"seqnames","chr")
   ref[,chr := as.character(chr)]
   ##TODO change to allow all contigs
   ref <- ref[chr %in% c(1:19,"X","Y")]
   
   setkey(ref,chr ,start,end)

   #read known_edit_file
   known_edit_tab <- fread(known_edit_file)
   setnames(known_edit_tab,c("chr","ref_pos"))
   known_edit_tab[,chr := as.character(chr)]
   setkey(known_edit_tab,chr,ref_pos)
   known_edit_tab_for_cluster_overlap <- known_edit_tab[,.(chr,cluster_start_pos = ref_pos,cluster_end_pos = ref_pos)]
   setkey(known_edit_tab_for_cluster_overlap)
   
   #read snp_file
   snp <- fread(snp_file, sep = "\t")
   setnames(snp,c("chr","ref_pos","pop_AF"))
   setkey(snp,chr,ref_pos)

   # #to test single chromosome  
   # x <- 1
   # chrom_ref = ref[chr == x,]
   # chrom_name = x
   # chrom_known_edit_tab = known_edit_tab[chr == x,]
   # chrom_snp_tab = snp[chr == x,]

   
   return("OK")
}
  



mismatch_tab <- fread("RNA_editing_analysis/VAL01-A2G002.mismatch_tab.tsv")
mismatch_tab <- mismatch_tab[ref == "A" & alt == "G" | ref == "T" & alt == "C" ]

window <- 50

mismatch_tab[,start_pos := ref_pos - window / 2 + 1] 
mismatch_tab[,end_pos := ref_pos + window / 2] 

# 
# intervals <- GRanges(
#   seqnames = mismatch_tab$chr,
#   ranges = IRanges(start = mismatch_tab$start_pos, end = mismatch_tab$end_pos)
# )

# Join overlapping intervals using reduce
# joined_intervals <- reduce(intervals)

pos_tab <- mismatch_tab[,.(chr,ref_pos)]
pos_tab <- unique(pos_tab)

pos_tab <- pos_tab[,.(pos = unlist(lapply(ref_pos,function(x) return(seq(x - window / 2 + 1,x + window / 2,1))))),by = chr]
pos_tab <- unique(pos_tab)

setorder(pos_tab)

fwrite(pos_tab,"test_loci.tsv",col.names = F,sep = "\t")

mismatch_tab

# samples <- gsub(".bam","",list.files("/mnt/ssd/ssd_1/workspace/vojta/val_editing/input_files/mapped/",pattern = ".bam$"))
# sample <- tail(samples,1)

setwd("/mnt/ssd/ssd_1/workspace/vojta/editing_nn_model")

res <- lapply(samples,function(sample){
  args <- character(5)
  args[1] <- paste0("/mnt/ssd/ssd_1/workspace/vojta/editing_nn_model/input_files/mapped/",sample,".bam") 
  args[2]
  args[3] <- "/mnt/nfs/shared/CFBioinformatics/references_backup/mus_musculus/GRCm38.p6-93/annot/GRCm38.p6-93.gtf"
  args[4] <- "/mnt/nfs/shared/CFBioinformatics/references_backup/mus_musculus/GRCm38.p6-93/other/known_editing_sites/GRCm38.p6-93.csv"
  args[5] <- "/mnt/nfs/shared/CFBioinformatics/references_backup/mus_musculus/GRCm38.p6-93/other/snp/GRCm38.p6-93.all_snp.tsv"
  args[6] <- paste0("/mnt/ssd/ssd_1/workspace/vojta/val_editing/RNA_editing_analysis/",sample)

  print(sample)
  # run_all(args)
})

# print(res)
 
#Run as Rscript
# args <- commandArgs(trailingOnly = T)
# run_all(args)
 






# 
# sam <- mclapply(unique(ref$chr),read_single_chr_bam_as_sam,bam_file = bam_file,mc.preschedule = F,mc.silent = T,mc.cleanup = T,mc.cores = 16)
# sam <- rbindlist(sam)
# 
# all_reads_count <- nrow(unique(sam, by = c("qname","read","mapq")))
# # fwrite(sam,file = paste0(output_file_prefix,".test"), sep = "\t")
# # sam <- fread(paste0(paste0("/mnt/ssd/ssd_1/workspace/vojta/A2I_edit_EXO10_KD/RNA_editing_analysis/",samples[24]),".test"), sep = "\t")
# 
# sam <- gene_annotation_and_multimapped_split(sam,ref[,.(chr,start,end,gene_id,strand)],iterations = 3,min_ratio_to_keep = 0.05)
# 
# mismatch_tab <- mclapply(unique(ref$chr),function(x) {
#    return(process_single_chromosome(chrom_sam = sam[chr == x,],
#                                     chrom_known_edit_tab = known_edit_tab[chr == x,],
#                                     chrom_snp_tab = snp[chr == x,]))
# },mc.preschedule = F,mc.silent = T,mc.cleanup = T,mc.cores = 16)
# 
# mismatch_tab <- rbindlist(mismatch_tab)
# 
# fwrite(mismatch_tab,file = paste0(output_file_prefix,".all_mismatch_tab.tsv"), sep = "\t")
# 
# # mismatch_tab <- fread("/mnt/ssd/ssd_1/workspace/vojta/A2I_edit_EXO10_KD/RNA_editing_analysis2/u87_adar_ctrl_rep1.all_mismatch_tab.tsv")
# 
# 
# #basic process + quality and population filter
# edit_signal_vars <- strsplit(strsplit(edit_setings,split = "\\|")[[1]][1],split = ";")[[1]]
# edit_signal_vars <- matrix(unlist(tstrsplit(edit_signal_vars,split = "")[c(2,4)]),nrow = 2)
# add_edit_signal_true_text <- paste0("mismatch_tab[",paste("(ref == '",edit_signal_vars[,1],"' & alt == '",edit_signal_vars[,2],"')",sep = "",collapse = " | "),",edit_signal := T]")
# 
# mismatch_tab[,edit_signal := F]
# eval(parse(text = add_edit_signal_true_text))
# mismatch_tab <- mismatch_tab[qual_prob > 96]
# mismatch_tab <- mismatch_tab[pop_AF < 0.01]
# mismatch_tab[,self_variant_count := .N,by = .(row_id,ref,alt)]
# known_editing_random_error_count <- nrow(known_edit_tab) * (1 - mean(mismatch_tab[edit_signal == T & known_editting == T,]$qual_prob)/100)
# 
# #edit regions
# var_dist_window <- 150
# editing_cluster_tab <- mclapply(unique(ref$chr),function(x) {
#    return(mismatch_tab[chr == x,single_var_tab_cluster_computation(.SD,var_dist_window),by = .(ref,alt)])
# },mc.preschedule = F,mc.silent = T,mc.cleanup = T,mc.cores = 16)
# 
# editing_cluster_tab <- rbindlist(editing_cluster_tab)
# setkey(editing_cluster_tab,chr,cluster_start_pos,cluster_end_pos)
# 
# known_edit_tab_overlap <- foverlaps(known_edit_tab_for_cluster_overlap,editing_cluster_tab,nomatch = NULL)
# known_edit_tab_overlap <- known_edit_tab_overlap[,.(known_editing_pos = .N),by = .(chr,cluster_start_pos,cluster_end_pos)]
# 
# editing_cluster_tab <- merge(editing_cluster_tab,known_edit_tab_overlap,by = c("chr","cluster_start_pos","cluster_end_pos"),all.x = T)
# editing_cluster_tab[is.na(known_editing_pos),known_editing_pos := 0]
# editing_cluster_tab <- merge(ref[,.(gene_id,gene_name,strand)],editing_cluster_tab,by = "gene_id",all.y = T)
# 
# fwrite(editing_cluster_tab,file = paste0(output_file_prefix,".editing_cluster_tab.tsv"), sep = "\t")

# test_res <- lapply(seq(from = 1,to = 20,by = 1),function(x) {
#    test <- editing_cluster_tab[var_count > 7 + 8 & position_count > 7 & mean_vars_per_read > 2]
#    return(test[,.(.N,sum(position_count)),by = .(ref,alt)])
# })
# 
# names(test_res) <- seq(from = 1,to = 20,by = 1)
# 
# lapply(test_res,function(x) sum(x[ref == "A" & alt == "G" | ref == "T" & alt == "C"]$N) / sum(x[ref == "G" & alt == "A" | ref == "C" & alt == "T"]$N))

# test[,coef4 := known_editing_var_pos / known_editing_pos]
# hist(test[ref == "A" & alt == "G"]$coef4,probability = T)
# test[,coef4 := known_editing_var_pos / known_editing_pos]
# test[,coef := cluster_size / position_count]
# test[,coef2 := position_count / var_count]
# hist(test[ref == "G" & alt == "A"]$mean_vars_per_read,probability = T)
# hist(test[ref == "A" & alt == "G"]$mean_vars_per_read,probability = T)
# hist(test[ref == "G" & alt == "A" & coef < 200]$coef2,probability = T)
# hist(test[ref == "A" & alt == "G" & coef < 200]$coef2,probability = T)
# hist(test[ref == "G" & alt == "A" & coef3 < 200]$coef3,probability = T)
#
# hist(test[ref == "A" & alt == "G"]$coef4,probability = T)
# hist(test[coef3 < 1000]$coef3,probability = T)
# test[,.(.N,sum(var_count)),by = .(ref,alt)]
#
# test[,.(.N,
#         mean(mean_vars_per_read),
#         mean(cluster_size / position_count),
#         mean(var_count - position_count)),by = .(ref,alt)]



# 
# mismatch_tab <- mismatch_tab[per_map_mismatch_count - self_variant_count <= 3]
# mismatch_tab <- mismatch_tab[self_variant_count >= 2 | known_editting == T]
# 
# fwrite(mismatch_tab,file = paste0(output_file_prefix,".mismatch_tab.tsv"), sep = "\t")
# 
# per_sample_res <- mismatch_tab[,.(all_reads_count = all_reads_count,
#                                   all_mapping_count = sum(sam$read_split_ratio),
#                                   all_vars = sum(read_split_ratio),
#                                   abs_edit_signal_level = sum(read_split_ratio * edit_signal),
#                                   known_editting = sum(read_split_ratio * known_editting),
#                                   edit_clusters = sum(editing_cluster_tab[var_count > 7 + 8 & position_count > 7 & mean_vars_per_read > 2]$edit_signal),
#                                   neg_clusters = sum(!editing_cluster_tab[var_count > 7 + 8 & position_count > 7 & mean_vars_per_read > 2]$edit_signal))]
# # per_sample_res <- per_sample_res[,.(all_vars = .N,abs_edit_signal_level = sum(edit_signal),known_editting = sum(known_editting)),by = .(sample)]
# per_sample_res[,abs_noise_level := all_vars - abs_edit_signal_level]
# per_sample_res[,corected_abs_edit_signal_level := abs_edit_signal_level - known_editing_random_error_count]
# per_sample_res[,edit_signal_noise_ratio := corected_abs_edit_signal_level / abs_noise_level]
# per_sample_res[,edit_clusters_read_ratio := edit_clusters / all_reads_count * 10^6]
# per_sample_res[,edit_signal_read_ratio := corected_abs_edit_signal_level / all_reads_count]
# per_sample_res[,edit_signal_mapping_ratio := corected_abs_edit_signal_level / all_mapping_count]
# per_sample_res[,log2_signal_noise := log2(edit_signal_noise_ratio)]
# per_sample_res[,editing_level := log2_signal_noise / 2.55 * 100]
# 
# fwrite(per_sample_res,file = paste0(output_file_prefix,".sample_editing.tsv"), sep = "\t")
# 
# per_gene_res <-  mismatch_tab[,.(all_vars = sum(read_split_ratio),abs_edit_signal_level = sum(read_split_ratio * edit_signal),known_editting = sum(read_split_ratio * known_editting)),by = gene_id]
# per_gene_res <- per_gene_res[all_vars > min_gene_vars]
# per_gene_res[,abs_noise_level := all_vars - abs_edit_signal_level]
# per_gene_res[,abs_edit_signal_level := abs_edit_signal_level + smooth_coef]
# per_gene_res[,abs_noise_level := abs_noise_level + smooth_coef]
# per_gene_res[,edit_signal_noise_ratio := abs_edit_signal_level / abs_noise_level]
# per_gene_res[,log2_signal_noise := log2(edit_signal_noise_ratio)]
# per_gene_res[,editing_level := log2_signal_noise / 2.55 * 100]
# per_gene_res[log2_signal_noise < 0,editing_level := 0]
# 
# 
# 
# 
# per_gene_res <- merge(ref[,.(gene_id,gene_name,strand)],per_gene_res,by = "gene_id")
# 
# fwrite(per_gene_res,file = paste0(output_file_prefix,".per_gene_editing.tsv"), sep = "\t")







 ######
 #####
 #####
 # number of A-G variants and all variants per read (add pair-end read mate information (not neccesary now))
 # mismatch_tab[ref == "A" & alt == "G", AG_variants := .N, by = mapping_id]
 # mismatch_tab[is.na(AG_variants),AG_variants := 0]
 # mismatch_tab[, all_variants := .N, by = mapping_id]
 
 #  mismatch_pos_MD <- lapply(str_match,function(x) cumsum(as.numeric(x[,2]) + 1) )
 #  orig_nuc <- lapply(str_match,function(x) x[,3] )
 #  rm(str_match)
 #  gc()
 #  num_of_mismatch <- sapply(mismatch_pos_MD,length)
 # # finding real positions of mismatches regarding cigar information
 #  tab_MD <- data.table(mp = rep(overlapped$mapping_id, times = num_of_mismatch), rel_pos = unlist(mismatch_pos_MD))
 #  if(length(tab_MD$mp)>0){
 #    tab_cigar <- data.table(mp = rep(overlapped$mapping_id, times = num_of_positions), num_SIM = unlist(num_of_SIM), positions = unlist(SIM_pos), SIM = unlist(S_or_I_or_M))
 #    joined <- merge(tab_MD, tab_cigar, by = c("mp"), all.x = TRUE)
 #    #filtering softClips and inserions as relative positions from MD are only affected by them.
 #    join_SI <- joined[SIM %like% "[SI]"]
 #    join_SI <- join_SI[,rel_pos := ifelse(.I[1], ifelse(rel_pos >= positions, rel_pos + num_SIM, rel_pos), ifelse(lag(rel_pos) >= positions, lag(rel_pos) + num_SIM, lag(rel_pos))), by = c("mp", "rel_pos")]
 #    join_M_other_reads <- joined[is.na(match(joined$mp, join_SI$mp)),]
 #    real_pos <- rbind(join_SI, join_M_other_reads)
 #    setkey(real_pos, mp)
 #    real_pos <- real_pos[,list(max.score=max(rel_pos)), by=c("mp", "rel_pos")]
 #    mismatch_tab <- data.table(mapping_id = real_pos$mp,real_pos = real_pos$max.score, ref = unlist(orig_nuc),
 #                           alt = stringi::stri_sub(rep(overlapped$read,times = num_of_mismatch),from = real_pos$max.score,length = 1),
 #                           varq = stringi::stri_sub(rep(overlapped$mapq,times = num_of_mismatch),from = real_pos$max.score,length = 1))
 #    mismatch_tab <- cbind(overlapped[mismatch_tab$mapping_id,.(chr,start,end,i.start,gene_id,strand,read_split_ratio)],mismatch_tab)
 #    mismatch_tab[,pos := i.start + real_pos - 1]
 #    rm(overlapped)
 #    gc()
 
 
 
 # nuc_change <- structure(c("A","C","G","T"),names = c("G","T","A","C"))
 # mismatch_tab[,trans_ref := ref]
 # mismatch_tab[,trans_alt := alt]
 # mismatch_tab[strand == "-",trans_ref := nuc_change[ref]]
 # mismatch_tab[strand == "-",trans_alt := nuc_change[alt]]
 
 
#Run as Rscript
# args <- commandArgs(trailingOnly = T)
# run_all(args)
#or
#
#result <- read_SS(oneFile.loc = "/mnt/ssd/ssd_3/temp/ailar/output.txt",
 #                   separateFile = FALSE, withMFE = TRUE)
# b) count per gene "editing factor"
#             1) compare A-G variants to to non A-G
#             2) compare A-G variants to all reads coverage
#             3) other improvments that improve difference between u87_adar_kd and u87_adar_control


  
 #keep only right strandness for aligned reads
 #overlapped <- overlapped[strand_flag == "+" & strand == "+" | strand_flag == "-" & strand == "-"]
 #overlapped[, strand_flag := NULL]
 
  
  
 # #fwrite(sam,gsub("bam$","chr1_test.sam",bam_file))
 # # #to test
 # 
 # sam <- fread("/mnt/ssd/ssd_1/workspace/vojta/A2I_edit_EXO10_KD/test_tab_ctrl.tsv")
 # setnames(sam,c("qname", "flag", "chr", "start", "cigar", "read", "mapq", "MD"))
 # sam <- sam[chr == "22"]
 # sam[,end := start + nchar(read)]
 # sam[,chr := as.character(chr)]
 # #  dt_AG <- fread(known_edit_file)  
 # #  setkey(dt_AG)
 
 
 
 #  fwrite(dt_AG[tab_seq$dna %in% c("A","T")],file = known_edit_file,row.names = F,col.names = F,sep = "\t")
 
 # sam <- fread("/mnt/ssd/ssd_1/snakemake/stage105_RNA_edit_pipeline_test.second/input_files/mapped/u87_adar_kd_rep2.chr1_test.sam")
 
 #adding strandness information from flag
 #check_strandness <- function(x){
 #  if (samFlags(x)[8] == FALSE){
 #   return("+")
 #  }else{
 #   return("-")
 # }
 #}
 #sam <- sam[,strand_flag := sapply(flag,function(x) check_strandness(x))]

 # #secondary_structure
 # tic("secondary_structure_input_preparation")
 # #Input_preparation
 # tab <- mismatch_tab[AG_variants > 0]
 # tab[, pos_min_fifty := pos - 50]
 # tab[, pos_plus_fifty := pos + 50]
 # 
 # #I tried pos+-hundred, but rnafold took around ten min by using only two cpus that wont be efficient. (I am wondering if we can use all number of cpus; parallel.cores == -1)
 # genome <- BSgenome.Hsapiens.UCSC.hg38
 # tab_ranges <- reduce(GRanges( paste0("chr",tab$chr),
 #                      IRanges(start=tab$pos_min_fifty, end = tab$pos_plus_fifty),
 #                      strand=tab$strand ))
 #                           
 # #dt_granges <- tab[, as.data.table(reduce(GRanges(paste0("chr", chr),IRanges(start=pos_min_fifty, end=pos_plus_fifty), strand=strand))), by = chr]
 # #tab_ranges <- GRanges(dt_granges$seqnames, IRanges(start=dt_granges$start, end=dt_granges$end), strand=dt_granges$strand)
 # sort(tab_ranges, ignore.strand=TRUE)
 # tab_seq <- Views(genome,tab_ranges)
 # tab_seq<- setorder(as.data.table(tab_seq))
 # tab_seq[, id := seq_along(dna)]
 # tab_seq[, id := paste0("seq", id)]
 # seq <- tab_seq$dna
 # names(seq) <- tab_seq$id
 # dna <- DNAStringSet(seq)
 # writeXStringSet(dna, fasta)
 # toc(log=T)
 # cat(tic.log()[[length(tic.log())]],"\n", file = log_file, append = TRUE, sep="" )
 # 
 # # secondary_structure_run
 # tic("secondary_structure_run")
 # seqs <- read.fasta(file = fasta)
 # SS.seq_2 <- run_RNAfold(seqs, RNAfold.path = "RNAfold", parallel.cores = 2)
 # toc(log=T)
 # cat(tic.log()[[length(tic.log())]],"\n", file = log_file, append = TRUE, sep="" )
 # 
 # tic("overlapping")
 # ss_seq <- as.data.table(t(SS.seq_2))
 # tab_seq <- cbind(tab_seq, ss_seq[, V2, V3])
 # tab_pos <- tab[, .(seqnames = chr, start = pos, end = pos, width = 1, strand, mapping_id)]
 # tab_seq[, seqnames := gsub("chr", "", seqnames)]
 # setkeyv(tab_seq, c("seqnames","start","end"))
 # overlapped_rna_fold <- foverlaps(tab_pos, tab_seq, by.x=c("seqnames","start","end"), by.y=c("seqnames","start","end"), nomatch = 0)
 # 
 # overlapped_rna_fold[, ss := stringi::stri_sub(V2,from = i.start - start + 1,length = 1)]
 # setnames(overlapped_rna_fold, "V3", "MFE")
 # mismatch_tab_edit <- merge(mismatch_tab, overlapped_rna_fold[, .(mapping_id, ss, MFE)], by = c("mapping_id"), all.x= TRUE)
 # toc(log=T)
 # cat(tic.log()[[length(tic.log())]],"\n", file = log_file, append = TRUE, sep="" )
 # write.table(mismatch_tab_edit, file = editing_sites, row.names = F, quote = F)
 #    }else{
 #      header_tab <- c("mapping_id", "chr", "pos", "start", "end",
 #                      "i.start", "gene_id", "strand", "read_split_ratio",
 #                      "real_pos", "ref", "alt", "P_accuracy", "known_editting",
 #                      "AG_variants", "all_variants", "pop_AF", "ss", "MFE")
 #      write(header_tab, editing_sites, ncolumns = 17, sep="\t")
 #      write("", fasta)
 #    }




# 
# 
# 
# 
# res <- mclapply(unique(ref$chr),function(x) {
#   return(process_single_chromosome(bam_file = bam_file,
#                                    chrom_ref = ref[chr == x,],
#                                    chrom_name = x,
#                                    chrom_known_edit_tab = known_edit_tab[chr == x,],
#                                    chrom_snp_tab = snp[chr == x,],
#                                    edit_variation = edit_variation,
#                                    var_dist_window = var_dist_window,
#                                    cluster_window_extension = cluster_window_extension))
# },mc.preschedule = F,mc.silent = T,mc.cleanup = T,mc.cores = 16)
# 
# edit_cluster = rbindlist(lapply(res,function(x) x$edit_cluster))
# mismatch_tab = rbindlist(lapply(res,function(x) x$mismatch_tab))
# 
# fwrite(edit_cluster,file = paste0(output_file_prefix,".editing_clusters.tsv"), sep = "\t")
# fwrite(mismatch_tab,file = paste0(output_file_prefix,".mismatch_tab.tsv"), sep = "\t")
# 
# # edit_cluster <- fread(paste0(output_file_prefix,".editing_clusters.tsv"))
# # gene_coverage_tab <- fread(paste0(output_file_prefix,".gene_coverage_tab.tsv"))
# 
# 
# 
# filtered_edit_cluster_tab <- edit_cluster[position_count > 3 & 
#                                             var_count > 4 & 
#                                             var_count / all_vars > 0.7 & 
#                                             position_count / all_vars_pos > 0.6 & 
#                                             (all_vars_extended / all_vars) / ((cluster_size + 2 * cluster_window_extension) / cluster_size) < 0.8 &
#                                             mean_median_bias < 2]
# 
# 
# fwrite(filtered_edit_cluster_tab,file = paste0(output_file_prefix,".filtered_editing_clusters.tsv"), sep = "\t")
# # 
# # edit_cluster[position_count > 10 & var_count > 20 & var_count / all_vars > 0.8 & read_basepairs / cluster_size > 6]
# # 
# # fitered_edit_cluster_tab <- edit_cluster[position_count > 2 & var_count > 4]
# # 
# # per_sample_res <- gene_coverage_tab[,.(all_basepair_count = sum(gene_coverage_tab$read_basepairs),
# #                                                                     all_mapping_count = sum(sam$read_split_ratio),
# #                                                                     all_vars = sum(read_split_ratio),
# #                                                                     abs_edit_signal_level = sum(read_split_ratio * edit_signal),
# #                                                                     known_editting = sum(read_split_ratio * known_editting),
# #                                                                     edit_clusters = sum(editing_cluster_tab[var_count > 7 + 8 & position_count > 7 & mean_vars_per_read > 2]$edit_signal),
# #                                                                     neg_clusters = sum(!editing_cluster_tab[var_count > 7 + 8 & position_count > 7 & mean_vars_per_read > 2]$edit_signal))]
# 
# toc(log=T)
# cat(tic.log()[[length(tic.log())]],"\n", file = paste0(output_file_prefix,".log"), append = TRUE, sep="" )

