library(data.table)
library(openxlsx)


script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(paste0(script_dir,"/.."))


old_code <- function(){
  cluster_window_extension = 20
  samples <- gsub(".bam","",list.files("input_files/mapped/",pattern = ".bam$"))
  
  # all_data_labels <- data.table(
  #   condition = gsub("_rep.*","",list.files("/mnt/ssd/ssd_1/snakemake/stage9998_OConnell_rna_editing_EXOSC10/input_files/per_sample_results/",pattern = ".editing_sites.all.tsv")),
  #   rep = gsub(".*_rep(.).*","rep_\\1",list.files("/mnt/ssd/ssd_1/snakemake/stage9998_OConnell_rna_editing_EXOSC10/input_files/per_sample_results/",pattern = ".editing_sites.all.tsv")),
  #   file = list.files("/mnt/ssd/ssd_1/snakemake/stage9998_OConnell_rna_editing_EXOSC10/input_files/per_sample_results/",pattern = ".editing_sites.all.tsv",full.names = T)
  # )
  # 
  # all_data_labels <- rbind(all_data_labels,data.table(
  #   condition = gsub("_rep.*","",list.files("/mnt/ssd/ssd_1/snakemake/stage9997_OConnell_rna_editing_DIS3/input_files/per_sample_results/",pattern = ".editing_sites.all.tsv")),
  #   rep = gsub(".*_rep(.).*","rep_\\1",list.files("/mnt/ssd/ssd_1/snakemake/stage9997_OConnell_rna_editing_DIS3/input_files/per_sample_results/",pattern = ".editing_sites.all.tsv")),
  #   file = list.files("/mnt/ssd/ssd_1/snakemake/stage9997_OConnell_rna_editing_DIS3/input_files/per_sample_results/",pattern = ".editing_sites.all.tsv",full.names = T)
  # ))
  # 
  # 
  # all_data_labels <- rbind(all_data_labels,data.table(
  #   condition = gsub("_rep.*","",list.files("/mnt/ssd/ssd_1/snakemake/stage105_RNA_edit_pipeline_test.second/input_files/per_sample_results/",pattern = ".editing_sites.all.tsv")),
  #   rep = gsub(".*_rep(.).*","rep_\\1",list.files("/mnt/ssd/ssd_1/snakemake/stage105_RNA_edit_pipeline_test.second/input_files/per_sample_results/",pattern = ".editing_sites.all.tsv")),
  #   file = list.files("/mnt/ssd/ssd_1/snakemake/stage105_RNA_edit_pipeline_test.second/input_files/per_sample_results/",pattern = ".editing_sites.all.tsv",full.names = T)
  # ))
  
  alu <- fread("input_files/Alu.bed")
  setnames(alu,c("chr","start","end","Alu_name","length","Alu_strand"))
  alu <- alu[grepl("alu",Alu_name,ignore.case = T)]
  alu[,chr := gsub("chr","",chr)]
  alu <- alu[chr %in% c(1:22,"X","Y")]
  alu <- alu[,.(chr,start,end,Alu_name,Alu_strand)]
  
  analysis_dir <- "RNA_editing_analysis"
  
  all_data_labels <- data.table(
    sample = samples,
    cluster_file = paste0(analysis_dir,"/",samples,".editing_clusters.tsv"),
    mismatch_tab_file = paste0(analysis_dir,"/",samples,".mismatch_tab.tsv")
  )
  
  sample_design <- fread("input_files/DE_experiment_design.tsv")
  sample_design[,batch_group := NULL]
  setnames(sample_design,"replicate","rep")
  setnames(sample_design,"sample_name","sample")
  
  all_data_labels <- merge(sample_design,all_data_labels,by = "sample")
  
  #variants_processing
  mismatch_tab <- lapply( all_data_labels$mismatch_tab_file,fread)
  names(mismatch_tab) <- all_data_labels$sample
  mismatch_tab <- rbindlist(mismatch_tab,use.names = T,idcol = "sample")
  
  per_sample_res <- mismatch_tab[,.(all_vars = .N),by = sample]
  remove(mismatch_tab)
  
  per_sample_res[,sample_normalization_factor := 1 / (all_vars / mean(all_vars))]
  
  #clusters
  cluster_res <- lapply(all_data_labels$cluster_file,fread)
  names(cluster_res) <- all_data_labels$sample
  cluster_res <- rbindlist(cluster_res,use.names = T,idcol = "sample")
  
  setorder(cluster_res,-var_count)
  
  
  #TODO - overlap with alu
  
  
  setkey(alu)
  setkey(cluster_res,chr,cluster_start_pos,cluster_end_pos)
  
  cluster_res <- foverlaps(cluster_res,alu,by.x = c("chr","cluster_start_pos","cluster_end_pos"),by.y = c("chr","start","end"))
  setnames(cluster_res,c("start","end"),c("Alu_start","Alu_end"))
  setcolorder(cluster_res,names(cluster_res)[c(1,6:20)])
  cluster_res <- unique(cluster_res,by = c("chr","sample","cluster_start_pos","cluster_end_pos"))
  
  #TODO - get coverage of all clusters (alus, genes)
  
  # nrow(cluster_res[is.na(Alu_name)]) / nrow(cluster_res)
  # nrow(filtered_edit_cluster_tab[is.na(Alu_name)]) / nrow(filtered_edit_cluster_tab)
  # nrow(filtered_edit_cluster_tab[!is.na(Alu_name)]) / nrow(cluster_res[!is.na(Alu_name)]) 
  # nrow(filtered_edit_cluster_tab) / nrow(cluster_res)
  # filtered_edit_cluster_tab[is.na(Alu_name)]
  
  
  #TODO - overlap between samples
  
  
  # filtered_edit_cluster_tab <- cluster_res[position_count > 3 & 
  #                                             var_count > 4 & 
  #                                             var_count / all_vars > 0.7 & 
  #                                             position_count / all_vars_pos > 0.6 & 
  #                                             (all_vars_extended / all_vars) / ((cluster_size + 2 * cluster_window_extension) / cluster_size) < 0.8 &
  #                                             mean_median_bias < 2]
  
  filtered_edit_cluster_tab <- cluster_res[position_count > 4 & 
                                             var_count > 12 & 
                                             var_count / all_vars > 0.7 & 
                                             position_count / all_vars_pos > 0.5 & 
                                             (all_vars_extended / all_vars) / ((cluster_size + 2 * cluster_window_extension) / cluster_size) < 0.8 &
                                             mean_median_bias < 2]
  
  
  sample_order <- setorder(filtered_edit_cluster_tab[,sum(cluster_end_pos - cluster_start_pos),by = sample],-V1)$sample
  
  
  add_sample <- sample_order[2]
  joined_tab <- filtered_edit_cluster_tab[sample == sample_order[1],.(chr,cluster_start_pos,cluster_end_pos,variant,gene_id,Alu_name,Alu_strand,var_count,sample)]
  joined_tab <- unique(joined_tab,by = c("chr","cluster_start_pos","cluster_end_pos","variant"))
  
  for(add_sample in sample_order[-1]){
    cut_tab <- joined_tab[,.(chr,cluster_start_pos,cluster_end_pos,variant,gene_id)]
    cut_tab <- unique(cut_tab)
    add_tab <- filtered_edit_cluster_tab[sample == add_sample,.(chr,cluster_start_pos,cluster_end_pos,variant,gene_id,Alu_name,Alu_strand,var_count,sample)]
    setkey(cut_tab,chr,cluster_start_pos,cluster_end_pos)
    setkey(add_tab,chr,cluster_start_pos,cluster_end_pos)
    
    add_tab <- foverlaps(add_tab,cut_tab)
    add_tab[is.na(cluster_start_pos),c("cluster_start_pos","cluster_end_pos","variant","gene_id") := .(i.cluster_start_pos,i.cluster_end_pos,i.variant,i.gene_id)]
    add_tab <- add_tab[,.(chr,cluster_start_pos,cluster_end_pos,variant,gene_id,Alu_name,Alu_strand,var_count,sample)] 
    add_tab <- unique(add_tab,by = c("chr","cluster_start_pos","cluster_end_pos","variant"))
    
    joined_tab <- rbind(joined_tab,add_tab)
  }
  
  
  joined_tab[,c("mean_vars","in_sample_count") := .(mean(var_count),.N),by = .(chr,cluster_start_pos,cluster_end_pos,variant,gene_id,Alu_name)]
  cast_joined_tab <- dcast.data.table(joined_tab,formula = chr + cluster_start_pos + cluster_end_pos + variant + gene_id + Alu_name + Alu_strand + mean_vars + in_sample_count ~ sample,value.var = "var_count",fill = 0)
  setorder(cast_joined_tab,-in_sample_count,-mean_vars)
  
  
  summary_results <- joined_tab[,.(clusters_raw = .N,vars_count_raw  = sum(var_count),non_unique_clusters_raw   = sum(in_sample_count > 1)),by = sample]
  
  # summary_results <- merge(summary_results,per_sample_res[,.(sample,sample_normalization_factor)],by = "sample")
  # summary_results[,clusters := clusters_raw * sample_normalization_factor]
  # summary_results[,vars_count := vars_count_raw * sample_normalization_factor]
  # summary_results[,non_unique_clusters := non_unique_clusters_raw * sample_normalization_factor]
  
  all_reads_tab <- fread("input_files/multiqc_star.txt") 
  all_reads_tab <- all_reads_tab[,.(sample = Sample,A_bp_count = (uniquely_mapped * avg_mapped_read_length)/4)]
  summary_results <- merge(summary_results,all_reads_tab,by = "sample")
  summary_results[,edited_per_1M_A_bp := vars_count_raw / A_bp_count * 10^6]
  setcolorder(summary_results,c("sample","edited_per_1M_A_bp"))
  summary_results <- merge(sample_design,summary_results,by = "sample")
  setorder(summary_results,condition,rep)
  
  
  write.xlsx(list(summary_results = summary_results,joined_clusters = cast_joined_tab),file ="val_editing_cluster_res.xlsx")
  
  
  
  
  
  
  
  test <- filtered_edit_cluster_tab[,.(clusters = .N,vars_count = sum(var_count),var_pos_count = sum(position_count)),by = sample]
  test <- merge(test,per_sample_res[,.(sample,sample_normalization_factor)],by = "sample")
  test[,clusters := clusters * sample_normalization_factor]
  test[,vars_count := vars_count * sample_normalization_factor]
  test[,var_pos_count := var_pos_count * sample_normalization_factor]
  
  strict_filtered_edit_cluster_tab <- cluster_res[position_count > 4 & 
                                                    var_count > 12 & 
                                                    var_count / all_vars > 0.7 & 
                                                    position_count / all_vars_pos > 0.5 & 
                                                    (all_vars_extended / all_vars) / ((cluster_size + 2 * cluster_window_extension) / cluster_size) < 0.8 &
                                                    mean_median_bias < 2]
  
  test2<- strict_filtered_edit_cluster_tab[,.(clusters = .N,vars_count = sum(var_count),var_pos_count = sum(position_count),cluster_size = mean(cluster_size),var_count = mean(var_count)),by = sample]
  test2 <- merge(test2,per_sample_res[,.(sample,sample_normalization_factor)],by = "sample")
  test2[,clusters := clusters * sample_normalization_factor]
  test2[,vars_count := vars_count * sample_normalization_factor]
  test2[,var_pos_count := var_pos_count * sample_normalization_factor]
  test2
  
  
  
  #per_sample
  sample_res <- lapply(all_data_labels$file,fread)
  names(sample_res) <- all_data_labels$sample
  sample_res <- rbindlist(sample_res,use.names = T,idcol = "sample")
  
  sample_res <- merge(all_data_labels[,.(condition,rep,sample)],sample_res,by = "sample")
  sample_res[,sample := NULL]
  setcolorder(sample_res,c("condition","rep","editing_level"))
  sample_res[,editing_level := round(editing_level,1)]
  
  vec <- (sample_res$all_vars / sample_res$all_mapping_count)
  data.table(all_data_labels$condition,all_data_labels$rep,vec / mean(vec))
  
  
  
  fwrite(sample_res,file = "/mnt/ssd/ssd_1/workspace/vojta/A2I_edit_EXO10_KD/RNA_editing_analysis/all_sample_editing.tsv",sep = "\t")
  
  # #per_gene
  # gene_res <- lapply( gsub(".sample_editing.tsv",".per_gene_editing.tsv",all_data_labels$file),fread)
  # names(gene_res) <- all_data_labels$sample
  # gene_res <- rbindlist(gene_res,use.names = T,idcol = "sample")
  # 
  # gene_res[,editing_level := round(editing_level,1)]
  # gene_res <- dcast.data.table(gene_res,gene_id + gene_name + strand ~ sample,value.var = "editing_level",fill = NA)
  # 
  # setorder(gene_res,gene_name)
  # 
  # fwrite(gene_res,file = "/mnt/ssd/ssd_1/workspace/vojta/A2I_edit_EXO10_KD/RNA_editing_analysis/per_gene_editing.tsv",sep = "\t")
  
  
  
  
  
  
  
  
  
  
  mismatch_tab <- mismatch_tab[self_variant_count >= 4 | (known_editting == T & self_variant_count >= 2)]
  
  
  per_sample_res <- mismatch_tab[,.(all_vars = .N,abs_edit_signal_level = sum(edit_signal),known_editting = sum(known_editting)),by = .(sample)]
  per_sample_res[,abs_noise_level := all_vars - abs_edit_signal_level]
  per_sample_res[,edit_signal_noise_ratio := abs_edit_signal_level / abs_noise_level]
  per_sample_res[,log2_signal_noise := log2(edit_signal_noise_ratio)]
  
  test_edit_tab[,.N,by = c("ref","alt")]
  # edit_tab[,.N,by = .(sample,ref,alt)]
  
  ######
  ######
  kd_test_tab <- edit_tab[sample %in% all_data_labels$sample[9:12]]
  kd_test_tab <- edit_tab
  kd_test_tab <- unique(kd_test_tab,by = c("sample","mapping_id","chr","pos"))
  kd_test_tab <- kd_test_tab[ref != "N" & alt != "N"]
  kd_test_tab <- kd_test_tab[ref != alt]
  kd_test_tab[,self_variant_count := .N,by = .(sample,mapping_id,ref,alt)]
  kd_test_tab[,all_variants_count := .N,by = .(sample,mapping_id)]
  
  kd_test_tab_bad <- kd_test_tab[all_variants_count - self_variant_count > 3]
  kd_test_tab <- kd_test_tab[all_variants_count - self_variant_count <= 3]
  
  kd_test_tab[,ref_type := 3]
  kd_test_tab[ref == "A" | ref == "T",ref_type := 2]
  kd_test_tab[,alt_type := 3]
  kd_test_tab[alt == "A" | alt == "T",alt_type := 2]
  
  
  control_set <- kd_test_tab[ref_type == alt_type]
  edit_set <- kd_test_tab[ref_type != alt_type]
  
  edit_set[((ref == "A" & alt == "G") | (ref == "T" & alt == "C")) & self_variant_count > 2][,.N,by = sample]
  control_set[self_variant_count > 2][,.N,by = sample]
  
  res <- merge(edit_set[((ref == "A" & alt == "G") | (ref == "T" & alt == "C")) & self_variant_count > 2][,.(edits = .N),by = sample],control_set[self_variant_count > 2][,.(norm = .N),by = sample],by = "sample")
  res[,edit_quotient := edits / norm]
  
  setnames(res,"edits","potential_A2I_edit_count")
  setnames(res,"norm","background_var_count")
  setnames(res,"edit_quotient","A2I_edit_index")
  setnames(res,"A2I_edit_index","relative_A2I_edit_index")
  
  res <- merge(all_data_labels,res,by = "sample")
  res[,sample := NULL]
  res[,file := NULL]
  
  res[,relative_A2I_edit_index := round(relative_A2I_edit_index,digits = 2)]
  
  write.xlsx(list(per_sample_editing = res),file = "first_look_edit_table.xlsx")
  
  
  # kd_test_tab <- kd_test_tab[all_variants < 4]
  kd_test_tab[,is_potential_edit := F]
  kd_test_tab[(ref == "A" & alt == "G") | (ref == "T" & alt == "C"),is_potential_edit := T]
  kd_test_tab[,is_AG := F]
  kd_test_tab[(ref == "A" & alt == "G") ,is_AG := T]
  kd_test_tab[,is_TC := F]
  kd_test_tab[(ref == "T" & alt == "C"),is_TC := T]
  
  kd_test_tab[,.N,.(ref,alt)]
  
  setorder(kd_test_tab,ref,alt)
  kd_test_tab[self_variant_count > 2][,.N,by = .(ref,alt,sample)]
  
  edit_set <- kd_test_tab[(ref == "A" & alt == "G") | (ref == "T" & alt == "C")]
  control_set <- kd_test_tab[(ref == "G" & alt == "A") | (ref == "C" & alt == "T")]
  
  control_set[,mean(self_variant_count),by = .(known_editting)]
  edit_set[,mean(self_variant_count),by = .(known_editting,)]
  
  edit_set[self_variant_count > 2][,.N,by = sample]
  control_set[self_variant_count > 2][,.N,by = sample]
  
  all_data_labels[,file := NULL]
  all_data_labels[,edit_count := sapply(res,nrow)]
  
  res2 <- lapply(res,function(tab) {
    tab <- tab[pop_AF < 0.02]
    tab <- tab[AG_variants > 1 & all_variants < 4 | known_editting == T]
    return(tab)
  })
  
  
  ###assign_cluster_var_size function
  # tab <- edit_set[sample == "u87_adar_ctrl__rep_1" & chr == "1" & ref == "A"]
  # pos <- tab$pos
  # 
  # assign_cluster_var_size <- function(pos,dist){
  #   sort_pos <- sort(pos)
  #   rle()
  #   diff(sort_pos)
  # }
  
  
  all_data_labels[,edit_count2 := sapply(res2,nrow)]
  
  all_data_labels
  
  table(tab$AG_variants)
  
}

