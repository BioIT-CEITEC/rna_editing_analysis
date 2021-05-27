library(stringi)
library(data.table)



run_all <- function(args){
  mpileup_file <- args[1]
  bam_file <- args[2]
  output_Rdata_file <- args[3]
  
  tab <- fread(mpileup_file)
  
  not_same_tab <- tab[V4 != nchar(V5),]
  not_same_tab[,V5 := gsub("\\^.","",V5)]
  not_same_tab[,V5 := gsub("\\$","",V5)]
  not_same_tab[,V5 := gsub("\\+[0-9]+[ACGTNacgtn]+","",V5)]
  not_same_tab[,V5 := gsub("-[0-9]+[ACGTNacgtn]+","",V5)]
  
  tab <- merge(tab,not_same_tab[,.(V1,V2,new_V5 = V5)],by = c("V1","V2"),all.x = T)
  tab[is.na(new_V5),new_V5 := V5]
  tab <- tab[V4 == nchar(new_V5),]
  
  tab[,mod_V5 := gsub("[\\<\\>]","",new_V5)]
  tab[,cover := nchar(mod_V5)]
  tab[,match := stringi::stri_count_regex(mod_V5,"[\\.\\,]")]
  tab[V3 == "A",other_mm := stringi::stri_count_regex(mod_V5,"[ATCatc]")]
  tab[V3 == "T",other_mm := stringi::stri_count_regex(mod_V5,"[ATGatg]")]
  per_pos_tab <- tab[,.(chrom = V1,pos = V2,ref = V3,cover,match,other_mm)]
  
  per_read_tab <- data.table(chrom = rep(tab$V1,tab$V4),pos = rep(tab$V2,tab$V4),read = unlist(strsplit(tab$V7,split = ",")),base = unlist(strsplit(tab$new_V5,split = "")))
  per_read_tab <- per_read_tab[base %in% c("G","g","N","n","C","c")]
  bam_file <-"/mnt/ssd/ssd_1/snakemake/stage105_RNA_edit_pipeline_test.second/input_files/mapped/u87_adar_kd_rep2.bam"
 
  sam <- fread(cmd=paste0("/mnt/ssd/ssd_3/install/anaconda3/bin/samtools view ",bam_file),fill = T,select = 1:20)
  setkey(sam,V1)
  mm_reads <- sam[unique(per_read_tab$read)]
  per_read_tab <- merge(per_read_tab,unique(setkey(mm_reads[,.(read = V1,read_mm = as.numeric(gsub(".*:","",V19)))]),by = "read"),by = "read")
  
  
  save(per_pos_tab,per_read_tab,file = output_Rdata_file)
  
}

# develop and test 
#args <- character(3)
#args[1] <- "/mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/input_files/results/u87_control_rep2.mpileup.tmp2.tsv" #mpileup_file
#args[2] <- "/mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/input_files/mapped/u87_control_rep2.sam" #bam_file
#args[3] <- "/mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/input_files/results/u87_control_rep2.edit_tab.Rdata" #output_file

 #run as Rscript
 args <- commandArgs(trailingOnly = T)
 run_all(args)














# library(GenomicRanges)
# library(BSgenome.Hsapiens.UCSC.hg38)
# 
# 
# 
# gtf_file <- "/mnt/ssd/ssd_3/references/homsap/GRCh38-p10/annot/GRCh38-p10.gtf"
# chr_list <- paste0("chr",c(1:22,"X","Y"))
# chr_list_NCBI <- c(1:22,"X","Y")
# 
# 
# genes_A_pos_from_GTF <- function(gtf_file){
#   gtf <- fread(gtf_file)  
#   gtf <- gtf[V3 == "gene"]
#   setnames(gtf,names(gtf)[1],"V1")
#   seqlevelsStyle(gtf$V1) <- "UCSC"
#   gtf <- gtf[V1 %in% chr_list,]
#   
#   
#   ranges <- GRanges( gtf$V1,IRanges(start=gtf$V4, end = gtf$V5),strand=gtf$V7 )
#   
#   seqs <- as.data.table(Views(Hsapiens,ranges))
#   seqs[,seqnames := as.character(seqnames)]
#   seqlevelsStyle(seqs$seqnames) <- "NCBI"
#   
#   pos <- stringi::stri_locate_all_fixed(seqs$dna,"A")
#   
#   res <- rbindlist(lapply(1:nrow(seqs),function(x){
#     data.table(chrom = seqs$seqnames[x],pos = pos[[x]][,1] + seqs$start[x] - 2)
#   }))
#   
#   res <- unique(res)
#   setkey(res)
#   fwrite(res,file = "/mnt/ssd/ssd_3/references/homsap/GRCh38-p10/intervals/RNA/all_As.pos",sep = "\t",col.names = F)
#   
#   res <- fread("/mnt/ssd/ssd_3/references/homsap/GRCh38-p10/intervals/RNA/all_As.pos")
# }
# 
# tab <- fread("results/079_post.mpileup_small.tsv")
# fwrite(tab[,.(V1,V2,V3,V4,V5)],file = "results/079_post.mpileup.small.tsv",sep = "\t",col.names = F)
# tab_filtered <- tab[grepl("[GgNn]",V5) & V3 == "A" | grepl("[CcNn]",V5) & V3 == "T",]
# 
# 
# 
# tab <- fread("unaligned_remapped/test.sam",fill=T)
# tab <- tab[V3 != "*",]
# tab <- tab[V3 %in% chr_list_NCBI]
# seqlevelsStyle(tab$V3) <- "UCSC"
# test_tab <- tab[grepl("^[0-9]+M",V6)]
# test_tab <- cbind(test_tab,as.data.table(matrix(unlist(as.binary(test_tab$V2,n = 9)),ncol = 9,byrow = T,dimnames = list(NULL,paste0("F",1:9)))))
# test_tab <- test_tab[F5 == F,]
# test_tab[,len := nchar(V10)]
# ranges <- GRanges( test_tab$V3,IRanges(start=test_tab$V4, end = test_tab$V4 + test_tab$len - 1),strand="+" )
# seqs <- as.data.table(Views(Hsapiens,ranges))
# test_tab <- cbind(test_tab,seqs$dna)

