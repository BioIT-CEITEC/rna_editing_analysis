#BiocManager::install("Rsamtools")
#install.packages("remotes")
#if(!require(SamSeq)) {remotes::install_github("jefferys/SamSeq");library(SamSeq)}
#BiocManager::install("IRanges")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#install.packages("seqinr")
if(!require(LncFinder)) {install.packages("LncFinder", repos = getCRANmirrors()$URL[1]); library("LncFinder")}
#install.packages("devtools")
#devtools::install_github("collectivemedia/tictoc")

library(IRanges)
library(data.table)
# library(Rsamtools)
library(rtracklayer)
library(stringr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(seqinr)
library(tictoc)

run_all <- function(args){
        bam_file <- args[1]
        gtf_file <- args[2]
        Human_AG_all_hg38 <- args[3]
        snp_file <- args[4]
        log_file <- args[5]
        fasta <- args[6]
        editing_sites <- args[7]

      sam <- fread(cmd = paste0("samtools view ", bam_file, " | cut -f 1-4,6,10,11,12"), col.names = c("qname", "flag", "chr", "start", "cigar", "read", "mapq", "MD"))
      sam <- sam[chr %in% c(1:21,"X","Y")]
   #take only first chromosome for test purposes 
  # sam <- sam[chr == "1"]
 
  #fwrite(sam,gsub("bam$","chr1_test.sam",bam_file))
  #to test
  #sam <- fread("/mnt/ssd/ssd_1/snakemake/stage105_RNA_edit_pipeline_test.second/input_files/mapped/u87_adar_kd_rep2.chr1_test.sam")
    sam[,end := start + nchar(read)]
    sam[,chr := as.character(chr)]
  #adding strandness information from flag
  #check_strandness <- function(x){
  #  if (samFlags(x)[8] == FALSE){
  #   return("+")
  #  }else{
 #   return("-")
  # }
 #}
  #sam <- sam[,strand_flag := sapply(flag,function(x) check_strandness(x))]
  #removing intergenic regions
   tic("removing_intergenic_regions")
   feat_type <- "gene"
   annotate_by<- c("gene_name","gene_id", "strand")  
   ref <- as.data.table(rtracklayer::import(gtf_file))[type == feat_type, c("seqnames","start","end",annotate_by), with=F]
   setkeyv(ref, c("seqnames","start","end"))
   overlapped <- foverlaps(sam, ref, by.x=c("chr","start","end"), by.y=c("seqnames","start","end"), nomatch = 0)
   rm(ref, sam)
   gc()
  #keep only right strandness for aligned reads
  #overlapped <- overlapped[strand_flag == "+" & strand == "+" | strand_flag == "-" & strand == "-"]
  #overlapped[, strand_flag := NULL]
   toc(log=T)
   cat(tic.log()[[length(tic.log())]],"\n", file = log_file, append = TRUE, sep="" )
   #log_open(log_file, append = TRUE)
   #log_print("removing_intergenic_regions :", exectime1)
   #log_close()
   tic("multimapped_reads_spliting_ratio")
   overlapped[,mapping_count := .N, by = .(qname,read,mapq)]
   overlapped[,normalized_reads_per_gene := .N /(end - start), by = gene_id]

# Calculating multimapped reads splitting ratio based on the gene expresions
   overlapped[, read_split_ratio := normalized_reads_per_gene/sum(normalized_reads_per_gene), by = .(qname,read,mapq)]
   toc(log=T)
   cat(tic.log()[[length(tic.log())]],"\n", file = log_file, append = TRUE, sep="" )
   #log_open(log_file, append = TRUE)
   #log_print("multimapped_reads_spliting_ratio :", exectime2)
   #log_close()
#getting cigar information
# I am removing hardClipping, introns and deleting information from both cigar and MD to get the right positions for insertions.
#we have asigned basepairs in MD as many as the numer for deletion in cigar. 
   tic("alternations")
   overlapped[, cigar_m := gsub("([0-9]+)([DHN])", "", cigar)]
   str_match_cigar <- str_match_all(overlapped$cigar_m, "([0-9]+)([SMI]+)")
   num_of_SIM <- lapply(str_match_cigar,function(x) as.numeric(x[,2]))
   SIM_pos <- lapply(str_match_cigar,function(x) cumsum(as.numeric(x[,2])) - as.numeric(x[,2]) +1)
   num_of_positions <- sapply(SIM_pos,length)
   S_or_I_or_M <- lapply(str_match_cigar,function(x) x[,3] )
   rm(str_match_cigar)
  #MD_information
   overlapped[, MD := gsub("MD:Z:", "", MD)]
   overlapped[,MD := ifelse(MD %like% "^[ACTG]", paste0('0', MD),MD)]
  #removing basepairs related to deletion from MD as mentioned before.
   overlapped[, MD := gsub("\\^[ACGTacgt]+","",MD)]
   overlapped[,mapping_id := seq_along(read)]
   setkey(overlapped,mapping_id)
   str_match <- str_match_all(overlapped$MD, "([0-9]+)([ACTGactg]+)")
   missmatch_pos_MD <- lapply(str_match,function(x) cumsum(as.numeric(x[,2]) + 1) )
   orig_nuc <- lapply(str_match,function(x) x[,3] )
   rm(str_match)
   gc()
   num_of_missmatch <- sapply(missmatch_pos_MD,length)
  # finding real positions of missmatches regarding cigar information
   tab_MD <- data.table(mp = rep(overlapped$mapping_id, times = num_of_missmatch), rel_pos = unlist(missmatch_pos_MD))
   tab_cigar <- data.table(mp = rep(overlapped$mapping_id, times = num_of_positions), num_SIM = unlist(num_of_SIM), positions = unlist(SIM_pos), SIM = unlist(S_or_I_or_M))
   join <- merge(tab_MD, tab_cigar, by = c("mp"), all.x = TRUE)
#filtering softClips and inserions as relative positions from MD are anly affected by them.
   join_SI <- join[SIM %like% "[SI]"]
   join_SI <- join_SI[,rel_pos := ifelse(.I[1], ifelse(rel_pos >= positions, rel_pos + num_SIM, rel_pos), ifelse(lag(rel_pos) >= positions, lag(rel_pos) + num_SIM, lag(rel_pos))), by = c("mp", "rel_pos")]
   join_M_other_reads <- join[is.na(match(join$mp, join_SI$mp)),]
   real_pos <- rbind(join_SI, join_M_other_reads)
   setkey(real_pos, mp)
   real_pos <- real_pos[,list(max.score=max(rel_pos)), by=c("mp", "rel_pos")]
   missmatch_tab <- data.table(mapping_id = real_pos$mp,real_pos = real_pos$max.score, ref = unlist(orig_nuc),
                            alt = stringi::stri_sub(rep(overlapped$read,times = num_of_missmatch),from = real_pos$max.score,length = 1),
                            varq = stringi::stri_sub(rep(overlapped$mapq,times = num_of_missmatch),from = real_pos$max.score,length = 1))
   missmatch_tab <- cbind(overlapped[missmatch_tab$mapping_id,.(chr,start,end,i.start,gene_id,strand,read_split_ratio)],missmatch_tab)
   missmatch_tab[,pos := i.start + real_pos - 1]
   rm(overlapped)
   gc()


   nuc_change <- structure(c("A","C","G","T"),names = c("G","T","A","C"))
  missmatch_tab[strand == "-",ref := nuc_change[ref]]
  missmatch_tab[strand == "-",alt := nuc_change[alt]]
  toc(log=T)
  cat(tic.log()[[length(tic.log())]],"\n", file = log_file, append = TRUE, sep="" )
  #log_open(log_file, append = TRUE)
  #log_print("alternations :", exectime3)
  #log_close()

# adding 'probability score' of real editing):
  tic("base_qualiy_score")
  missmatch_tab[, Q := as.numeric(sapply(varq, function(x) as.integer(charToRaw(as.character(x)))-33))]
  missmatch_tab[, P_error := 10^-(Q/10)]
  missmatch_tab[, P_accuracy := (1-P_error) * 100]
  missmatch_tab[,c("varq","Q","P_error") := NULL]
  toc(log=T)
  cat(tic.log()[[length(tic.log())]],"\n", file = log_file, append = TRUE, sep="" )
  #log_open(log_file, append = TRUE)
  #log_print("base_quality_score :", exectime4)
  #log_close()

#Finding known editing sites
#Genome lifting
#Human_AG_all_hg19 <- "/mnt/ssd/ssd_3/references/homsap/RNA_editing_AG_all/Human_AG_all_hg19_v2.txt"   (#check if this txt file is zero or one based)
#dt_AG <- fread(Human_AG_all_hg19)

#dt_AG <- dt_AG[, .(chromosome, position)]
#dt_AG[, start := position - 1]
#dt_AG[,end := position]
#dt_AG[,position := NULL]
#dt_AG$end <- paste0("-",dt_AG$end)
#dt_AG$start_end = paste0(dt_AG$start,dt_AG$end)
#dt_AG[, c("start", "end") := NULL]
#write.table(dt_AG,"/mnt/ssd/ssd_3/temp/ailar/AG.bed", sep = ":", row.names = F, quote = F, col.names = F)
  tic("known_editing_sites")

  dt_AG <- fread(Human_AG_all_hg38)
  dt_AG <- dt_AG[, .(chr = V1, pos = V3)]
  dt_AG$chr <- gsub("chr","",dt_AG$chr)
  setkey(dt_AG, chr)
  dt_AG[, known_editting := TRUE]
  missmatch_tab <- merge(x = missmatch_tab, y = dt_AG, by = c("chr", "pos"), all.x = TRUE)
  missmatch_tab[is.na(known_editting),known_editting := F]
  toc(log=T)
  cat(tic.log()[[length(tic.log())]],"\n", file = log_file, append = TRUE, sep="" )
  #log_open(log_file, append = TRUE)
  #log_print("known_editing_sites :", exectime5)
  #log_close()
  rm(dt_AG)
  gc()
# number of A-G variants and all variants per read (add pair-end read mate information (not neccesary now))
  missmatch_tab[ref == "A" & alt == "G", AG_variants := .N, by = mapping_id]
  missmatch_tab[is.na(AG_variants),AG_variants := 0]
  missmatch_tab[, all_variants := .N, by = mapping_id]

# check snips
  tic("snp")
  snp <- fread(snp_file, sep = "\t")
  snp <- snp[,c("V1", "V3", "V9")]
  setnames(snp, c("V1", "V3", "V9"), c("chr", "pos", "pop_AF"))
  snp[,pop_AF := as.numeric(str_extract(pop_AF, "(?<=;AF=)[0-9]+(.[0-9]+)?"))]
  snp$chr <- as.character(snp$chr)
  missmatch_tab <- merge(x = missmatch_tab, y = snp, by = c("chr", "pos"), all.x = TRUE)
  missmatch_tab[is.na(pop_AF),pop_AF := 0]
  toc(log=T)
  cat(tic.log()[[length(tic.log())]],"\n", file = log_file, append = TRUE, sep="" )
  #log_open(log_file, append = TRUE)
  #log_print("snp :", exectime6)
  #log_close()
  rm(snp)
  gc()

#secondary_structure
  tic("secondary_structure_input_preparation")
#Input_preparation
  tab <- missmatch_tab[AG_variants > 0]
  tab[, pos_min_fifty := pos - 50]
  tab[, pos_plus_fifty := pos + 50]

#I tried pos+-hundred, but rnafold took around ten min by using only two cpus that wont be efficient. (I am wondering if we can use all number of cpus; parallel.cores == -1)
   genome=BSgenome.Hsapiens.UCSC.hg38
   tab_ranges <- reduce(GRanges( paste0("chr",tab$chr),
                           IRanges(start=tab$pos_min_fifty, end = tab$pos_plus_fifty),
                           strand=tab$strand ))
                                
   #dt_granges <- tab[, as.data.table(reduce(GRanges(paste0("chr", chr),IRanges(start=pos_min_fifty, end=pos_plus_fifty), strand=strand))), by = chr]
   

   #tab_ranges <- GRanges(dt_granges$seqnames, IRanges(start=dt_granges$start, end=dt_granges$end), strand=dt_granges$strand)
  sort(tab_ranges, ignore.strand=TRUE)
  tab_seq <- Views(genome,tab_ranges)
  tab_seq<- setorder(as.data.table(tab_seq))
  tab_seq[, id := seq_along(dna)]
  tab_seq[, id := paste0("seq", id)]
  seq = tab_seq$dna
  names(seq) = tab_seq$id
  dna = DNAStringSet(seq)
  writeXStringSet(dna, fasta)
  toc(log=T)
  cat(tic.log()[[length(tic.log())]],"\n", file = log_file, append = TRUE, sep="" )
  #log_open(log_file, append = TRUE)
  #log_print("ss_input_preparaion :", exectime7)
  #log_close()
# secondary_structure_run
  tic("secondary_structure_run")
  seqs <- read.fasta(file = fasta)
  SS.seq_2 <- run_RNAfold(seqs, RNAfold.path = "RNAfold", parallel.cores = 2)
  toc(log=T)
  cat(tic.log()[[length(tic.log())]],"\n", file = log_file, append = TRUE, sep="" )
  #log_open(log_file, append = TRUE)
  #log_print("ss_run :", exectime8)
  #log_close()
  tic("overlapping")
  ss_seq <- as.data.table(t(SS.seq_2))
  tab_seq <- cbind(tab_seq, ss_seq[, V2, V3])
  tab.gr <- GRanges(IRanges(start=tab$pos, end = tab$pos), seqnames = tab$chr, strand = tab$strand)
  tab_pos <- as.data.table(tab.gr)
  tab_pos <- cbind(tab_pos, tab[, mapping_id])
  tab_seq[, seqnames := gsub("chr", "", seqnames)]
  setkeyv(tab_seq, c("seqnames","start","end"))
  overlapped_rna_fold <- foverlaps(tab_pos, tab_seq, by.x=c("seqnames","start","end"), by.y=c("seqnames","start","end"), nomatch = 0)

  overlapped_rna_fold[, ss := stringi::stri_sub(V2,from = i.start - start + 1,length = 1)]
  setnames(overlapped_rna_fold, c("i.V2", "V3"), c("mapping_id", "MFE"))
  missmatch_tab_edit <- merge(missmatch_tab, overlapped_rna_fold[, .(mapping_id, ss, MFE)], by = c("mapping_id"), all.x= TRUE)
  toc(log=T)
  cat(tic.log()[[length(tic.log())]],"\n", file = log_file, append = TRUE, sep="" )
  write.table(missmatch_tab_edit, file = editing_sites, row.names = F, quote = F)
}

    #args <- character(7)
    #args[1] <-"/mnt/ssd/ssd_1/snakemake/stage105_RNA_edit_pipeline_test.second/input_files/mapped/u87_adar_kd_rep2_adar_kd_rep2.bam"
    #args[2] <- "/mnt/ssd/ssd_3/references/homsap/GRCh38-p10/annot/GRCh38-p10.gtf"
    #args[3] <- "/mnt/ssd/ssd_3/references/homsap/GRCh38-p10/other/known_editing_sites/Human_AG_all_hg38_v2.csv"
    #args[4] <- "/mnt/ssd/ssd_3/references/homsap/GRCh38-p10/other/snp/GRCh38-p10.snp.bed"
    #args[5] <- "/mnt/ssd/ssd_1/snakemake/stage105_RNA_edit_pipeline_test.second/RNA_editing_analysis/sample_logs/u87_adar_kd_rep2/get_editing_sites"
    #args[6] <- "/mnt/ssd/ssd_3/temp/ailar/seq.fa"
    #args[7] <- "/mnt/ssd/ssd_3/temp/ailar/editing_sites.tsv"

#Run as Rscript
args <- commandArgs(trailingOnly = T)
run_all(args)
#or
#
#result <- read_SS(oneFile.loc = "/mnt/ssd/ssd_3/temp/ailar/output.txt",
 #                   separateFile = FALSE, withMFE = TRUE)
# b) count per gene "editing factor"
#             1) compare A-G variants to to non A-G
#             2) compare A-G variants to all reads coverage
#             3) other improvments that improve difference between u87_adar_kd and u87_adar_control


  
  
  
  



