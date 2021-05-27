library(data.table)
library(stringi)

run_all <- function(args){
  mpileup_file <- args[1]
  output_pos_file <- args[2]
  
  tab <- fread(mpileup_file)
  
  tab[,actual_cov := V4 - stri_count_fixed(V5,"<") - stri_count_fixed(V5,">")]
  fwrite(tab_filtered[,.(V1,V2,actual_cov)],file = gsub("\\.filtered_pos\\.",".all_pos_cov.",output_pos_file),sep = "\t",col.names = F)
  
  tab_filtered <- tab[V4 > 1,]
  tab_filtered <- tab[grepl("[GgNn]",V5) & V3 == "A" | grepl("[CcNn]",V5) & V3 == "T",]
  
  fwrite(tab_filtered[,.(V1,V2)],file = output_pos_file,sep = "\t",col.names = F)
  
}

# develop and test 
# args <- character(2)
# args[1] <- "/mnt/ssd/ssd_3/references/homsap/GRCh38-p10/annot/vep/" #mpileup_file
# args[2] <- "muta_taster" #output_pos_file


# run as Rscript
args <- commandArgs(trailingOnly = T)
run_all(args)
