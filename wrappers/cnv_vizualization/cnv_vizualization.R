suppressMessages(library(data.table))
suppressMessages(library(ggplot2))

run_all <- function(args){
  cnv_tabs_rdata <- args[1]
  out_filename_pdf <- args[2]
  
  load(cnv_tabs_rdata) 
  
  ggplot(data=cnv_intervals[sample %in% unique(cnv_intervals$sample)[1:16] & chr == "1"], aes(x=snp_pos, y=copy_estimate, color=sample)) +
    geom_smooth(method = "loess")

  ggplot(data=cnv_intervals[sample %in% unique(cnv_intervals$sample)[1:8] & chr == "1"], aes(x=snp_pos, y=copy_estimate_rollmean_50, group=sample,color = sample)) +
    geom_line()
  
  ggplot(data=cnv_intervals[sample %in% unique(cnv_intervals$sample)[1:2] & chr == "1"], aes(x=snp_pos, y=het_sum_50, group=sample,color = sample)) +
    geom_line()
  
  dev.off()
}


# develop and test
args <- "/mnt/ssd/ssd_1/snakemake/stage172_test_CNV_analysis/CNV_analysis/results/CNV_tabs.Rdata"
args <- c(args,"/mnt/ssd/ssd_3/references/homsap/GRCh38-p10/other/cnv_intervals/all_samples_cnv.pdf")


#run as Rscript
# 
# script_dir <- dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))
# args <- commandArgs(trailingOnly = T)
# run_all(args)



