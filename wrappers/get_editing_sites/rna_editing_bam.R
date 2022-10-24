library(data.table,verbose = F)
library(rtracklayer,verbose = F)
library(stringr,verbose = F)
library(parallel,verbose = F)
options(scipen=999) #disable scientific notation (because we are creating BED file)

read_single_chr_bam_as_sam <- function(bam_file, chr_name){

  sam <- fread(cmd = paste0("samtools view ",bam_file," ",chr_name," | awk -v OFS='\t' '$0 !~ /MD:Z:[0-9]+\t/{match($0,/MD:Z:[^\t]+/); print $1, $2, $3, $4, $6, $10, $11, substr($0,RSTART,RLENGTH)}'"), col.names = c("qname", "flag", "chr", "start", "cigar", "read", "mapq", "MD"))

  sam[,end := start + nchar(read)]
  sam[,chr := as.character(chr)]

  return(sam)
}


gene_annotation_and_multimapped_split <- function(sam,chr_ref){

  sam[,row_id := seq_along(chr)]
  sam <- foverlaps(sam,chr_ref)
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

  # mismatch_tab[,alt := stringr::str_sub( sam[J(mismatch_tab$row_id)]$read,read_pos,read_pos)]
  # mismatch_tab[,base_quality := stringr::str_sub( sam[J(mismatch_tab$row_id)]$mapq,read_pos,read_pos)]

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
}


process_single_chromosome <- function(bam_file,chr_ref,chr_name){

  sam <- read_single_chr_bam_as_sam(bam_file,chr_name)
  sam <- gene_annotation_and_multimapped_split(sam,chr_ref[,.(chr,start,end,gene_id,strand)])
  cigar_tab <- create_cigar_table(sam)
  mismatch_tab <- create_mismatch_table(cigar_tab,sam)

  setkey(mismatch_tab,chr,ref_pos)

  return(mismatch_tab)
}


mismatch_coverage_count <- function(bam_file,mismatch_tab,output_file_prefix){
  mismatch_tab <- mismatch_tab[,.N,.(chr,ref_pos,ref,alt)][N > 1]

  pos_tab <- unique(mismatch_tab[,.(chr,ref_pos)])
  pos_tab[,ref_pos_end := ref_pos]
  pos_tab[,ref_pos := ref_pos-1]
  # fwrite(pos_tab,file = paste0(output_file_prefix,".positions_tmp.bed"), sep = "\t", col.names = F)
  fwrite(pos_tab,file = output_file_prefix, sep = "\t", col.names = F)

  # cov_tab <- fread(cmd = paste0("bedtools coverage -a ",output_file_prefix,".positions_tmp.bed -b ",bam_file," -counts"), col.names = c("chr", "start", "end", "cov"))
  cov_tab <- fread(cmd = paste0("bedtools coverage -a ",output_file_prefix," -b ",bam_file," -counts"), col.names = c("chr", "start", "end", "cov"))
  ## to test
  # cov_tab <- fread("bara_RNA_editing/Adar_1.2L.bed", col.names = c("chr", "start", "end", "cov"))

  mismatch_cov_tab <- merge(mismatch_tab,cov_tab,by.x=c("chr","ref_pos"),by.y=c("chr","end"))
  mismatch_cov_tab$start <- NULL

  # file.remove(paste0(output_file_prefix,".positions_tmp.bed"))

  return(mismatch_cov_tab)
}


run_all <- function(args){

  bam_file <- args[1]
  gtf_file <- args[2]
  output_file_prefix <- args[3]
  output_file <- args[4]

  # output_file_prefix <- gsub(".mismatch_cov_tab.tsv"," ",output_file)


  #read gtf_file
  feat_type <- "gene"
  annotate_by<- c("gene_name","gene_id", "strand")
  ref <- as.data.table(rtracklayer::import(gtf_file, feature.type = feat_type))[, c("seqnames","start","end",annotate_by), with=F]
  setnames(ref,"seqnames","chr")
  ref[,chr := as.character(chr)]
  ref <- ref[chr %in% c("2L","2R","3L","3R","4","X","Y")]
  setkey(ref,chr ,start,end)


  res <- mclapply(unique(ref$chr),function(x) {
    return(process_single_chromosome(bam_file = bam_file,
                                     chr_ref = ref[chr == x,],
                                     chr_name = x))
  },mc.preschedule = F,mc.silent = T,mc.cleanup = T,mc.cores = 16)
  mismatch_tab <- rbindlist(res)

  mismatch_cov_tab <- mismatch_coverage_count(bam_file,mismatch_tab,output_file_prefix)
  fwrite(mismatch_cov_tab,file = output_file, sep = "\t")
}



args <- commandArgs(trailingOnly = T)
# to test
# args <- c("/mnt/nfs/shared/RNA_and_Immunity/sequencing_results/projects/M6A_KD_flies_editing/RNA_alignment/mapped/Adar_1.bam",
#           "/mnt/ssd/ssd_3/references/drosophila_melanogaster/BDGP6-99/annot/BDGP6-99.gtf",
#           "Adar_1",
#           "/mnt/ssd/ssd_1/workspace/katka/bara_RNA_editing/Adar_1.mismatch_cov_tab.tsv")
run_all(args)

