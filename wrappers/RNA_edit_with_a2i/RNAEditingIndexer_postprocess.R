library(data.table)

test <- fread("/mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/a2i_test/RNA_editing_analysis_with_a2i/EditingIndex.csv")

test[,c("Sample",names(test)[grep("EditingIndex",names(test))]),with = F]

