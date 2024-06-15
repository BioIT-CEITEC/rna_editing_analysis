library(data.table)
library(tidyverse)
library(JACUSA2helper)
library(plyranges)
library(magrittr)

potential_edit_sites <- function(args){

    input_dir <- args[1]
    output_dir <- args[2]

    files <- list.files(path = input_dir, pattern = ".txt$", full.names = TRUE)

    for (file in files) {
        sample <- gsub(".txt", "", basename(file))
        
        raw_data_tab <- read_result(file)
        
        raw_data_tab$bc <- lapply_cond(raw_data_tab$bases, function(b) { Reduce("+", b) } ) %>% Reduce("+", .) %>% base_count()
        
        raw_data_tab <- raw_data_tab %>%
            dplyr::filter(score >= 2) %>%
            dplyr::filter(All(cov$cond1 >= 10)) %>%
            dplyr::filter(bc <= 2) %>%
            dplyr::filter(robust(bases))
            
            rna_bases_data <- Reduce("+", raw_data_tab$bases$cond1)
            ref2rna_data <- base_sub(rna_bases_data, raw_data_tab$ref)
            results_data_tab <- as.data.table(raw_data_tab)
            results_data_tab <- cbind(results_data_tab, ref2rna_data)
        
        results_data_tab[, sample := sample]
        results_data_tab[,c("info", "filter", "width", "bc", "FALSE.") := NULL]

        fwrite(results_data_tab, file = paste0(output_dir, "/", sample, "_potential_edit_sites.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    }
}