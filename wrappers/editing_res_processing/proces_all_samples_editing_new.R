library(data.table)
library(openxlsx)


script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(paste0(script_dir,"/../.."))


read_gtf <- function(gtf_file = "/mnt/share/share/710000-CEITEC/713000-cmm/713016-bioit/resources/references/drosophila_melanogaster/BDGP6/annot/r111/BDGP6.gtf"){
  feat_type <- "exon"
  annotate_by<- c("gene_name","gene_id", "strand","exon_number")
  ref <- as.data.table(rtracklayer::import(gtf_file, feature.type = feat_type))
  ref[,exon_version := as.integer(exon_version)]
  ref[,transcript_width := sum(width),by = transcript_id]
  setorder(ref,transcript_support_level,-transcript_version,-transcript_width,seqnames,start,exon_number,na.last = T)
  ref_tr_ids <- ref[,.(select_transcript_id = transcript_id[1]),by = gene_id]
  select_tr_ref <- ref[transcript_id %in% ref_tr_ids$select_transcript_id ]
  ref <- select_tr_ref[, c("seqnames","start","end",annotate_by), with=F]
  setnames(ref,"seqnames","chr")
  ref[,chr := as.character(chr)]
  ##TODO change to allow all contigs
  ref <- ref[chr %in% c("3R","3L","2R","2L","4","X","Y")]
  
  setkey(ref,chr ,start,end,strand)
  return(ref)
}

create_know_editing_sites_set <- function(gtf_tab){
  
  mismatch_tab <- fread("")
  read_gtf
  
  known_sites_editing_tab <- lapply(list.files("input_files/editing_analysis/known_editing_sites/",full.names = T),fread)
  names(known_sites_editing_tab) <-  gsub(".known_editing_sites_AF.tsv","",list.files("input_files/editing_analysis/known_editing_sites/"))
  known_sites_editing_tab <- rbindlist(known_sites_editing_tab,use.names = T,idcol = "sample")
  # known_sites_editing_tab <- known_sites_editing_tab[Good_depth > 0]
  # known_sites_editing_tab[,strand := "+"]
  # known_sites_editing_tab[Count_T + Count_C > Count_G + Count_A,strand := "-"]
  
  known_sites_editing_tab <- known_sites_editing_tab[,.(sample = sample,chr = `#CHR`,start = POS,end = POS,edit_count = Count_C + Count_G,cov = Good_depth)]
  
  setkey(known_sites_editing_tab,chr,start,end)
  setkey(gtf_tab,chr,start,end)
  known_sites_editing_tab <- foverlaps(known_sites_editing_tab,gtf_tab[,.(chr,start,end,gene_name,gene_id,exon_number)],by.y = c("chr","start","end"))
  
  
  known_sites_editing_tab[,pos := i.start]
  known_sites_editing_tab[,c("i.start","i.end","start","end") := NULL]
  setcolorder(known_sites_editing_tab,c("sample","chr","pos"))
  
  return(known_sites_editing_tab)
  # res <- known_sites_editing_tab[,.(cov = sum(Good_depth),edit_count = sum(Count_G) + sum(Count_C)),by = sample]
  # res[,ratio := edit_count / cov * 100]
} 


process_annot_tab <- function(annot_file = "editing_analysis/potential_edit_sites.annotated.tsv"){
  annot_tab <- fread(annot_file,sep = "\t",header = T,skip = "#Uploaded_variation",verbose = F,showProgress = F)
  
  setnames(annot_tab,"#Uploaded_variation","var_name")
  
  annot_tab_extra_parse <- annot_tab$Extra
  annot_tab_extra_parse <- strsplit(annot_tab_extra_parse,";")
  annot_tab_extra_names_order <- order(unlist(lapply(annot_tab_extra_parse,seq_along)))
  annot_extra_names <- unlist(lapply(annot_tab_extra_parse,function(x) gsub("(.*)=.*","\\1",x)))
  annot_extra_value <- unlist(lapply(annot_tab_extra_parse,function(x) gsub(".*=(.*)","\\1",x)))
  annot_extra_index <- sapply(annot_tab_extra_parse,length)
  annot_extra_index <- rep(seq_along(annot_extra_index),annot_extra_index)
  annot_tab_extra <- data.table(index = annot_extra_index,names = annot_extra_names,value = annot_extra_value)
  annot_tab_extra <- dcast.data.table(annot_tab_extra,formula = index ~ names,fill = NA,value.var = "value")
  
  annot_tab <- cbind(annot_tab,annot_tab_extra)
  remove(annot_tab_extra)
  remove(annot_tab_extra_parse)
  remove(annot_extra_names)
  remove(annot_extra_value)
  remove(annot_extra_index)
  remove(annot_tab_extra_names_order)
  
  annot_tab[,Extra := NULL]
  annot_tab[,index := NULL]
  
  annot_tab <- unique(annot_tab,by = c("var_name","Location","Allele","Gene","Feature","STRAND"))
  
  annot_tab[,c("chrom","pos","substitution") := as.data.table(tstrsplit(var_name,"_"))]
  annot_tab[,c("reference","alternative") := as.data.table(tstrsplit(substitution,"/"))]
  annot_tab[,pos := as.integer(pos)]
  
  annot_tab <- annot_tab[!is.na(STRAND)]
  annot_tab <- annot_tab[(STRAND == -1 & reference == "T") | (STRAND == 1 & reference == "A")]
  annot_tab[,int_dist := as.integer(DISTANCE)]
  annot_tab[is.na(int_dist),int_dist := 0]
  annot_tab[,int_can := 0]
  annot_tab[!is.na(CANONICAL),int_can := 1]
  setorder(annot_tab,-int_can,int_dist)
  annot_tab <- unique(annot_tab,by = c("var_name","BIOTYPE"))
  
  return(annot_tab)
}

know_editing_sites_analysis <- function(){
  
  annot_tab <- process_annot_tab()
  
  cast_annot_tab <- copy(annot_tab)
  cast_annot_tab[STRAND == 1,annot_STRAND := "+"]
  cast_annot_tab[STRAND == -1,annot_STRAND := "-"]
  cast_annot_tab[,annot := paste(SYMBOL,Gene,annot_STRAND,sep = " ")]
  cast_annot_tab <- dcast.data.table(data = cast_annot_tab,formula = chrom + pos ~ BIOTYPE,value.var = "annot",fill = "-")
  setcolorder(cast_annot_tab,c("chrom","pos","protein_coding","transposable_element","pseudogene","ncRNA"))
  
  editing_tab <- lapply(list.files("editing_analysis/potential_edit_site_AFs/",full.names = T),fread)
  names(editing_tab) <-  gsub(".editing_sites_AF.tsv","",list.files("editing_analysis/potential_edit_site_AFs/"))
  editing_tab <- rbindlist(editing_tab,use.names = T,idcol = "sample")
  setnames(editing_tab,c("#CHR","POS","Good_depth"),c("chrom","pos","cov"))
  editing_tab <- merge(editing_tab,unique(annot_tab[,.(chrom,pos,alternative)]),by = c("chrom","pos"))
  setorder(editing_tab,-cov)
  editing_tab <- unique(editing_tab,by = c("sample","chrom","pos"))
  editing_tab[alternative == "C",edit_count := Count_C]
  editing_tab[alternative == "G",edit_count := Count_G]
  editing_tab[,edit_ratio := edit_count / cov]
  
  
  editing_tab[,mean_cov := mean(cov),by = .(chrom,pos)]
  editing_tab <- editing_tab[mean_cov > 10,]

  
  editing_tab[,cond := gsub("_.$","",sample)]
  editing_tab[,rep := gsub("^.*_","",sample)]
  editing_tab <- editing_tab[,.(cond,rep,chrom,pos,edit_count,cov,edit_ratio)]
  
  
  annot_editing_tab <- merge(editing_tab,cast_annot_tab,by = c("chrom","pos"))
  
  
  per_sample_res <- annot_editing_tab[,.(mean_edit_ratio = mean(edit_ratio,na.rm = T),
                                   mean_edit_ratio_in_protCode = mean(edit_ratio[protein_coding != "-"],na.rm = T),
                                   mean_edit_ratio_in_TE = mean(edit_ratio[transposable_element != "-"],na.rm = T),
                                   mean_edit_ratio_in_pseudo = mean(edit_ratio[pseudogene != "-"],na.rm = T)),by = .(cond,rep)]
  
  setorder(per_sample_res,cond,rep)
  
  
  annot_editing_tab[,mean_gene_cov := round(mean(cov),1),by = .(protein_coding,transposable_element)]
  annot_editing_tab[,mean_gene_edit := round(mean(edit_ratio,na.rm = T),3),by = .(protein_coding,transposable_element)]
  per_gene_res <- annot_editing_tab[protein_coding != "-",.(transposable_element = transposable_element[1],mean_edit_ratio = round(mean(edit_ratio,na.rm = T),2)),by = .(cond,rep,protein_coding,mean_gene_cov,mean_gene_edit)]
  per_gene_res <- dcast.data.table(per_gene_res,protein_coding + transposable_element + mean_gene_cov + mean_gene_edit ~ cond + rep,value.var = "mean_edit_ratio",fill = NA)
  setorder(per_gene_res,-mean_gene_edit,-mean_gene_cov)
  
  editing_tab[,edit_ratio_round := round(edit_ratio,3)]
  all_edit_sites_tab <- dcast.data.table(editing_tab,chrom + pos ~ cond + rep,value.var = "edit_ratio_round",fill = NA)
  all_edit_sites_tab <- merge(cast_annot_tab,all_edit_sites_tab,by = c("chrom","pos"))
  
  fwrite(annot_editing_tab,file = "editing_analysis/full_annot_editing_tab.tsv")
  
  write.xlsx(per_sample_res,file = "editing_analysis/per_sample_editing.xlsx")
  write.xlsx(per_gene_res,file = "editing_analysis/per_coding_gene_editing.xlsx")
  write.xlsx(all_edit_sites_tab,file = "editing_analysis/all_editing_sites_tab.xlsx")
  
  
  setkey(known_sites_editing_tab,chr,start,end)
  setkey(gtf_tab,chr,start,end)
  known_sites_editing_tab <- foverlaps(known_sites_editing_tab,gtf_tab[,.(chr,start,end,gene_name,gene_id,exon_number)],by.y = c("chr","start","end"))
  
  
  known_sites_editing_tab[,pos := i.start]
  known_sites_editing_tab[,c("i.start","i.end","start","end") := NULL]
  setcolorder(known_sites_editing_tab,c("sample","chr","pos"))
  
  return(known_sites_editing_tab)
  # res <- known_sites_editing_tab[,.(cov = sum(Good_depth),edit_count = sum(Count_G) + sum(Count_C)),by = sample]
  # res[,ratio := edit_count / cov * 100]
} 

load_deseq2_tab <- function(DE_dir_name = "DE_RSEM",count_type = "normCounts"){
  
  DE_dir_path <- paste0("input_files/",DE_dir_name)
  de_tab_dirs <- list.files(DE_dir_path,pattern = "_vs_")
  
  tab <- NULL
  for(selected_de_tab_dir in de_tab_dirs){
    tmp_tab <- fread(paste0(DE_dir_path,"/",selected_de_tab_dir,"/DESeq2.tsv"))
    if(is.null(tab)){
      tmp_tab <- tmp_tab[,c("Ensembl_Id","pvalue",grep(count_type,names(tmp_tab),value = T)),with = F]
    } else {
      tmp_tab <- tmp_tab[,c("Ensembl_Id","pvalue"),with = F]
    }
    setnames(tmp_tab,"pvalue",paste0(selected_de_tab_dir,"_pval"))
    if(is.null(tab)){
      tab <- copy(tmp_tab)
    } else {
      tab <- merge(tab,tmp_tab,by ="Ensembl_Id")
    }
  }
  
  setcolorder(tab,c("Ensembl_Id",paste0(de_tab_dirs,"_pval")))
  
  return(tab)
}


run_all <- function(){
  DE_dir_name = "DE_RSEM"
  cluster_window_extension = 20
  samples <- gsub(".known_editing_sites_AF.tsv$","",list.files("input_files/editing_analysis/known_editing_sites/",pattern = ".known_editing_sites_AF.tsv$"))
  
  gtf_tab <- read_gtf()
  known_sites_editing_tab <- know_editing_sites_analysis(gtf_tab)
  expression_tab <- load_deseq2_tab()
  
  DE_dir_path <- paste0("input_files/",DE_dir_name)
  de_tab_dirs <- list.files(DE_dir_path,pattern = "_vs_")
  
  expression_melt_tab <- melt.data.table(expression_tab,id.vars = c("Ensembl_Id",paste0(de_tab_dirs,"_pval")),variable.name = "sample",value.name = "norm_read_count")
  expression_melt_tab[,sample := gsub("_normCounts","",sample)]
  expression_melt_tab[,sample := gsub("\\.","-",sample)]
  setnames(expression_melt_tab,"Ensembl_Id","gene_id")
  known_sites_editing_tab <- merge(known_sites_editing_tab,expression_melt_tab,by = c("sample","gene_id"))
  
  sample_design <- fread("input_files/DE_experiment_design.tsv")
  sample_design[,batch_group := NULL]
  sample_design[,condition2 := gsub("..$","",sample_name)]
  setnames(sample_design,"replicate","rep")
  setnames(sample_design,"sample_name","sample")
  
  known_sites_editing_tab <- merge(sample_design,known_sites_editing_tab,by = "sample")
  known_sites_editing_tab[,edit_ratio_perc_txt := paste0(round(edit_count / cov * 100,1)," (",edit_count,"/",cov,")"),]
  known_sites_editing_tab[,edit_ratio := edit_count / cov]
  known_sites_editing_tab[is.na(edit_ratio),edit_ratio := 0]
  known_sites_editing_tab[,sum_edit_counts := sum(edit_count),by = .(chr,pos)]
  known_sites_editing_tab[,sum_cov := sum(cov),by = .(chr,pos)]
  known_sites_editing_tab[,sample_desc := paste(condition2,condition,rep,sep = " ")]
  known_sites_editing_tab <- unique(known_sites_editing_tab,by = c("chr","pos", "gene_name","sample_desc"))
  known_sites_editing_tab_cast <- dcast.data.table(known_sites_editing_tab,chr + pos + gene_name + exon_number + sum_cov + sum_edit_counts ~ sample_desc,value.var = "edit_ratio_perc_txt")
  known_sites_editing_tab_cast <- known_sites_editing_tab_cast[sum_edit_counts > 0]
  
  setorder(known_sites_editing_tab_cast,chr,pos)
  known_sites_editing_tab_cast_filtered <- known_sites_editing_tab_cast[sum_edit_counts > 100 ]
  
  
  # write.xlsx(known_sites_editing_tab_cast,file = "val_editing_sites_table.xlsx")
  # write.xlsx(known_sites_editing_tab_cast_filtered,file = "val_editing_sites_table_filtered.xlsx")
  # 
  # summary_sample <- known_sites_editing_tab_filtered[,.(edit_count = sum(edit_count),cov = sum(cov),editing_ratio = round(sum(edit_count) / sum(cov),3),over_50perc_edit_sites = sum(edit_ratio > 0.5)),by = .(sample,condition,rep)]
  # summary_ADAR_WT <- known_sites_editing_tab_filtered[,.(edit_count = sum(edit_count),cov = sum(cov),editing_ratio = round(sum(edit_count) / sum(cov),3),over_50perc_edit_sites = sum(edit_ratio > 0.5) / length(unique(sample))),by = .(condition)]
  # summary_PER_PL <- known_sites_editing_tab_filtered[,.(edit_count = sum(edit_count),cov = sum(cov),editing_ratio = round(sum(edit_count) / sum(cov),3),over_50perc_edit_sites = sum(edit_ratio > 0.5) / length(unique(sample))),by = .(condition2)]
  # summary_ADAR_WT_PER_PL <- known_sites_editing_tab_filtered[,.(edit_count = sum(edit_count),cov = sum(cov),editing_ratio = round(sum(edit_count) / sum(cov),3),over_50perc_edit_sites = sum(edit_ratio > 0.5) / length(unique(sample))),by = .(condition,condition2)]
  # 
  # write.xlsx(list(summary_sample,summary_ADAR_WT,summary_PER_PL,summary_ADAR_WT_PER_PL),file = "editing_summary.xlsx")
  
  # write.xlsx(known_sites_editing_tab_cast_filtered,file = "val_editing_sites_table_filtered.xlsx")
  
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
  
  filtered_edit_cluster_tab <- cluster_res[position_count > 2 & 
                                             var_count > 4 & 
                                             var_count / all_vars > 0.7 & 
                                             position_count / all_vars_pos > 0.6 & 
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
  
  
  # joined_tab <- joined_tab[,.(var_count = sum(var_count)),by = .(chr,cluster_start_pos,cluster_end_pos,variant,gene_id, Alu_name, Alu_strand,sample)]
  joined_tab[,c("mean_vars","in_sample_count") := .(mean(var_count),.N),by = .(chr,cluster_start_pos,cluster_end_pos,variant,gene_id,Alu_name)]
  joined_tab <- merge(joined_tab,unique(gtf_tab[,.(gene_id,gene_name)],by = "gene_id"),by = "gene_id",all.x = T)
  joined_tab[gene_id == "ncRNA",gene_name := "NA"]
  joined_tab[gene_id == "ncRNA",gene_id := "NA"]
  cast_joined_tab <- dcast.data.table(joined_tab,formula = chr + cluster_start_pos + cluster_end_pos + variant + gene_id + gene_name + Alu_name + Alu_strand + mean_vars + in_sample_count ~ sample,value.var = "var_count",fill = 0)
  setorder(cast_joined_tab,-in_sample_count,-mean_vars)
  
  cast_joined_tab_cond <- copy(cast_joined_tab)
  setnames(cast_joined_tab_cond,sample_design$sample,paste(sample_design$condition,sample_design$rep,sep = "_"))
  
  by_gene_results <- joined_tab[gene_id != "NA",.(cluster_edits_count_raw  = sum(var_count)),by = .(sample,gene_id,gene_name)]
  
  known_sites_editing_tab_by_gene <- known_sites_editing_tab[gene_id != "NA" & sum_edit_counts > 50,.(known_edit_sites = .N,cov = sum(cov),editing_ratio = sum(edit_count)/sum(cov)),by = .(sample,gene_id,gene_name)]
  by_gene_results <- merge(known_sites_editing_tab_by_gene,by_gene_results,by = c("sample","gene_id","gene_name"),all = T)
  
  
  expression_melt_tab <- melt.data.table(expression_tab,id.vars = c("Ensembl_Id","ADAR_vs_WT_pval"),variable.name = "sample",value.name = "norm_read_count")
  expression_melt_tab[,sample := gsub("_normCounts","",sample)]
  expression_melt_tab[,sample := gsub("\\.","-",sample)]
  setnames(expression_melt_tab,"Ensembl_Id","gene_id")
  by_gene_results <- merge(expression_melt_tab,by_gene_results,by = c("sample","gene_id"))
  
  raw_read_gene_tab <- load_deseq2_tab(count_type = "rawCounts")
  raw_read_gene_tab[,grep("_pval",names(raw_read_gene_tab),value = T) := NULL]
  raw_read_gene_tab <- melt.data.table(raw_read_gene_tab,id.vars = c("Ensembl_Id"),variable.name = "sample",value.name = "raw_read_count")
  raw_read_gene_tab[,sample := gsub("_rawCounts","",sample)]
  raw_read_gene_tab[,sample := gsub("\\.","-",sample)]
  setnames(raw_read_gene_tab,"Ensembl_Id","gene_id")
  by_gene_results <- merge(by_gene_results,raw_read_gene_tab,by = c("sample","gene_id"))
  
  
  by_gene_results <- merge(sample_design,by_gene_results,by = "sample")
  
  setorder(by_gene_results,ADAR_vs_WT_pval,condition,rep,na.last = T)
  by_gene_results <- by_gene_results[!is.na(ADAR_vs_WT_pval)]
  by_gene_results[,cluster_edits_per_10000_reads := cluster_edits_count_raw / (raw_read_count + 10) * 10000]
  setcolorder(by_gene_results,c("gene_id","gene_name","ADAR_vs_WT_pval","sample","condition","rep","norm_read_count","raw_read_count","editing_ratio"))
  
  by_sample_results <- by_gene_results[,.(editing_ratio_perc = sum(editing_ratio * cov,na.rm = T) / sum(cov,na.rm = T) * 100,cluster_edits_per_10000_reads = sum(cluster_edits_count_raw,na.rm = T) / sum(raw_read_count) * 10000),by = .(sample,condition,rep)]
  by_sample_results[,editing_ratio_perc := round(editing_ratio_perc,2)]
  by_sample_results[,cluster_edits_per_10000_reads := round(cluster_edits_per_10000_reads,2)]
  
  by_gene_results[is.na(editing_ratio), editing_ratio := NA]
  by_gene_results[is.na(known_edit_sites), known_edit_sites := 0]
  by_gene_results[,editing_ratio := round(editing_ratio*100,2)]
  by_gene_results[,cluster_edits_per_10000_reads := round(cluster_edits_per_10000_reads,2)]
  setnames(by_gene_results,"editing_ratio","editing_ratio_perc")
  
  # summary_results <- merge(summary_results,per_sample_res[,.(sample,sample_normalization_factor)],by = "sample")
  # summary_results[,clusters := clusters_raw * sample_normalization_factor]
  # summary_results[,vars_count := vars_count_raw * sample_normalization_factor]
  # summary_results[,non_unique_clusters := non_unique_clusters_raw * sample_normalization_factor]
  
  # all_reads_tab <- fread("input_files/multiqc_star.txt") 
  # all_reads_tab <- all_reads_tab[,.(sample = Sample,A_bp_count = (uniquely_mapped * avg_mapped_read_length)/4)]
  # summary_results <- merge(summary_results,all_reads_tab,by = "sample")
  # summary_results[,edited_per_1M_A_bp := vars_count_raw / A_bp_count * 10^6]
  # setcolorder(summary_results,c("sample","edited_per_1M_A_bp"))
  # summary_results <- merge(sample_design,summary_results,by = "sample")
  # setorder(summary_results,condition,rep)
  
  
  write.xlsx(list(by_sample = by_sample_results,by_gene = by_gene_results),file ="val_editing_new_res.xlsx")
  write.xlsx(list(condition_names = cast_joined_tab_cond,sample_names = cast_joined_tab),file ="val_editing_cluster_res2.xlsx")
  
}



test_tdg <- function(){
  
  analysis_dir <- "RNA_editing_analysis"
  samples <- gsub(".known_editing_sites_AF.tsv$","",list.files("input_files/editing_analysis/known_editing_sites/",pattern = ".known_editing_sites_AF.tsv$"))
  
  mismatch_tab_file = paste0(analysis_dir,"/",samples,".mismatch_tab.tsv")
  mismatch_tab <- lapply(mismatch_tab_file,fread)
  names(mismatch_tab) <- samples
  mismatch_tab <- rbindlist(mismatch_tab,use.names = T,idcol = "sample")
  
  tdg_missmatches <- mismatch_tab[fwd_gene_id == "ENSMUSG00000034674"]
  tdg_missmatches_A_G <- tdg_missmatches[ref == "A" & alt == "G"]
  counts <- tdg_missmatches_A_G[,.N, by = .(sample,ref_pos)]
  counts[,sum := sum(N),by = ref_pos]
  counts_cast <- dcast.data.table(counts,ref_pos + sum ~ sample,value.var = "N",fill = 0)
  setorder(counts_cast,-sum)
  counts_cast_cond <- copy(counts_cast)
  setnames(counts_cast_cond,sample_design$sample,paste(sample_design$condition,sample_design$rep,sep = "_"))
  write.xlsx(list(condition_names = counts_cast_cond,sample_names = counts_cast),file ="TDG_A_to_G_var_counts_per_position.xlsx")
  
  tdg_missmatches <- mismatch_tab[rev_gene_id == "ENSMUSG00000031328"]
  tdg_missmatches_A_G <- tdg_missmatches[ref == "A" & alt == "G"]
  counts <- tdg_missmatches_A_G[,.N, by = .(sample,ref_pos)]
  counts[,sum := sum(N),by = ref_pos]
  counts_cast <- dcast.data.table(counts,ref_pos + sum ~ sample,value.var = "N",fill = 0)
  setorder(counts_cast,-sum)
  counts_cast_cond <- copy(counts_cast)
  setnames(counts_cast_cond,sample_design$sample,paste(sample_design$condition,sample_design$rep,sep = "_"))
  write.xlsx(list(condition_names = counts_cast_cond,sample_names = counts_cast),file ="TDG_A_to_G_var_counts_per_position.xlsx")
  
  sample_design <- fread("input_files/DE_experiment_design.tsv")
  sample_design[,batch_group := NULL]
  sample_design[,condition2 := gsub("..$","",sample_name)]
  setnames(sample_design,"replicate","rep")
  setnames(sample_design,"sample_name","sample")
  
  mismatch_tab_test <- merge(sample_design,mismatch_tab,by ="sample")
  
  res_fwd <- mismatch_tab_test[ref == "A" & alt == "G" ,.N,by = .(condition,fwd_gene_id)]
  
  # test <- mismatch_tab_test[fwd_gene_id != "ncRNA" & ref == "A" & alt == "G",.N,by = .(fwd_gene_id)]
  
  res_fwd[,sum_N := sum(N),by = .(fwd_gene_id)]
  res_fwd <- res_fwd[sum_N > 400,N[2]/N[1],by = .(fwd_gene_id)]
  res_fwd <- res_fwd[!is.na(V1)]
  setorder(res_fwd,V1,na.last = T)
  filt_res_fwd <- res_fwd[log2(V1) < -1.5]
  
  res_rev <- mismatch_tab_test[ref == "T" & alt == "C" ,.N,by = .(condition,rev_gene_id)]
  
  # test <- mismatch_tab_test[fwd_gene_id != "ncRNA" & ref == "A" & alt == "G",.N,by = .(fwd_gene_id)]
  
  res_rev[,sum_N := sum(N),by = .(rev_gene_id)]
  res_rev <- res_rev[sum_N > 400,N[2]/N[1],by = .(rev_gene_id)]
  res_rev <- res_rev[!is.na(V1)]
  setorder(res_rev,V1,na.last = T)
  filt_res_rev <- res_rev[log2(V1) < -1.5]
  
  
  
  mismatch_tab_filt_genes_fwd <- mismatch_tab_test[fwd_gene_id %in% filt_res_fwd$fwd_gene_id,.(condition,chr,ref_pos,gene = fwd_gene_id,ref,alt,stran = "FWD")]
  mismatch_tab_filt_genes_rev<- mismatch_tab_test[rev_gene_id %in% filt_res_rev$rev_gene_id,.(condition,chr,ref_pos,gene = rev_gene_id,ref,alt,stran = "REV")]
  
  mismatch_tab_filt_genes <- rbind(mismatch_tab_filt_genes_fwd,mismatch_tab_filt_genes_rev)
  fwrite(mismatch_tab_filt_genes,"mismatch_tab_wired_genes.tsv")
  
  
  
  filt_genes_counts <- mismatch_tab_filt_genes[,.(count = .N),by = .(condition,chr,gene,ref_pos,ref,alt, stran)]
  filt_genes_counts <- dcast.data.table(filt_genes_counts,chr + gene + ref_pos + ref + alt + stran ~ condition,value.var = "count",fill = 0 )
  filt_genes_counts[,sum_count := ADAR + WT]
  
  threshold = 50
  
  filt_genes_counts <- filt_genes_counts[sum_count > threshold]
  filt_genes_counts[,l2_ratio_AD_WT := log2((WT+1)/(ADAR+1))]
  setorder(filt_genes_counts,l2_ratio_AD_WT)
  
  threshold = -2
  filt_changed_pos <- filt_genes_counts[l2_ratio_AD_WT < threshold]
  
  setorder(filt_changed_pos,gene,ref_pos,ref,alt)
  # filt_changed_pos[,test := .N,by = ref_pos]
  # filt_changed_pos[test > 1]
  

  

  
  # hist(log2(res2$V1),40)
  
}

# 
# test <- filtered_edit_cluster_tab[,.(clusters = .N,vars_count = sum(var_count),var_pos_count = sum(position_count)),by = sample]
# test <- merge(test,per_sample_res[,.(sample,sample_normalization_factor)],by = "sample")
# test[,clusters := clusters * sample_normalization_factor]
# test[,vars_count := vars_count * sample_normalization_factor]
# test[,var_pos_count := var_pos_count * sample_normalization_factor]
# 
# strict_filtered_edit_cluster_tab <- cluster_res[position_count > 4 & 
#                                            var_count > 12 & 
#                                            var_count / all_vars > 0.7 & 
#                                            position_count / all_vars_pos > 0.5 & 
#                                            (all_vars_extended / all_vars) / ((cluster_size + 2 * cluster_window_extension) / cluster_size) < 0.8 &
#                                            mean_median_bias < 2]
# 
# test2<- strict_filtered_edit_cluster_tab[,.(clusters = .N,vars_count = sum(var_count),var_pos_count = sum(position_count),cluster_size = mean(cluster_size),var_count = mean(var_count)),by = sample]
# test2 <- merge(test2,per_sample_res[,.(sample,sample_normalization_factor)],by = "sample")
# test2[,clusters := clusters * sample_normalization_factor]
# test2[,vars_count := vars_count * sample_normalization_factor]
# test2[,var_pos_count := var_pos_count * sample_normalization_factor]
# test2
# 
# 
# 
# #per_sample
# sample_res <- lapply(all_data_labels$file,fread)
# names(sample_res) <- all_data_labels$sample
# sample_res <- rbindlist(sample_res,use.names = T,idcol = "sample")
# 
# sample_res <- merge(all_data_labels[,.(condition,rep,sample)],sample_res,by = "sample")
# sample_res[,sample := NULL]
# setcolorder(sample_res,c("condition","rep","editing_level"))
# sample_res[,editing_level := round(editing_level,1)]
# 
# vec <- (sample_res$all_vars / sample_res$all_mapping_count)
# data.table(all_data_labels$condition,all_data_labels$rep,vec / mean(vec))
# 
# 
# 
# fwrite(sample_res,file = "/mnt/ssd/ssd_1/workspace/vojta/A2I_edit_EXO10_KD/RNA_editing_analysis/all_sample_editing.tsv",sep = "\t")
# 
# # #per_gene
# # gene_res <- lapply( gsub(".sample_editing.tsv",".per_gene_editing.tsv",all_data_labels$file),fread)
# # names(gene_res) <- all_data_labels$sample
# # gene_res <- rbindlist(gene_res,use.names = T,idcol = "sample")
# # 
# # gene_res[,editing_level := round(editing_level,1)]
# # gene_res <- dcast.data.table(gene_res,gene_id + gene_name + strand ~ sample,value.var = "editing_level",fill = NA)
# # 
# # setorder(gene_res,gene_name)
# # 
# # fwrite(gene_res,file = "/mnt/ssd/ssd_1/workspace/vojta/A2I_edit_EXO10_KD/RNA_editing_analysis/per_gene_editing.tsv",sep = "\t")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# mismatch_tab <- mismatch_tab[self_variant_count >= 4 | (known_editting == T & self_variant_count >= 2)]
# 
# 
# per_sample_res <- mismatch_tab[,.(all_vars = .N,abs_edit_signal_level = sum(edit_signal),known_editting = sum(known_editting)),by = .(sample)]
# per_sample_res[,abs_noise_level := all_vars - abs_edit_signal_level]
# per_sample_res[,edit_signal_noise_ratio := abs_edit_signal_level / abs_noise_level]
# per_sample_res[,log2_signal_noise := log2(edit_signal_noise_ratio)]
# 
# test_edit_tab[,.N,by = c("ref","alt")]
# # edit_tab[,.N,by = .(sample,ref,alt)]
# 
# ######
# ######
# kd_test_tab <- edit_tab[sample %in% all_data_labels$sample[9:12]]
# kd_test_tab <- edit_tab
# kd_test_tab <- unique(kd_test_tab,by = c("sample","mapping_id","chr","pos"))
# kd_test_tab <- kd_test_tab[ref != "N" & alt != "N"]
# kd_test_tab <- kd_test_tab[ref != alt]
# kd_test_tab[,self_variant_count := .N,by = .(sample,mapping_id,ref,alt)]
# kd_test_tab[,all_variants_count := .N,by = .(sample,mapping_id)]
# 
# kd_test_tab_bad <- kd_test_tab[all_variants_count - self_variant_count > 3]
# kd_test_tab <- kd_test_tab[all_variants_count - self_variant_count <= 3]
# 
# kd_test_tab[,ref_type := 3]
# kd_test_tab[ref == "A" | ref == "T",ref_type := 2]
# kd_test_tab[,alt_type := 3]
# kd_test_tab[alt == "A" | alt == "T",alt_type := 2]
# 
# 
# control_set <- kd_test_tab[ref_type == alt_type]
# edit_set <- kd_test_tab[ref_type != alt_type]
# 
# edit_set[((ref == "A" & alt == "G") | (ref == "T" & alt == "C")) & self_variant_count > 2][,.N,by = sample]
# control_set[self_variant_count > 2][,.N,by = sample]
# 
# res <- merge(edit_set[((ref == "A" & alt == "G") | (ref == "T" & alt == "C")) & self_variant_count > 2][,.(edits = .N),by = sample],control_set[self_variant_count > 2][,.(norm = .N),by = sample],by = "sample")
# res[,edit_quotient := edits / norm]
# 
# setnames(res,"edits","potential_A2I_edit_count")
# setnames(res,"norm","background_var_count")
# setnames(res,"edit_quotient","A2I_edit_index")
# setnames(res,"A2I_edit_index","relative_A2I_edit_index")
# 
# res <- merge(all_data_labels,res,by = "sample")
# res[,sample := NULL]
# res[,file := NULL]
# 
# res[,relative_A2I_edit_index := round(relative_A2I_edit_index,digits = 2)]
# 
# write.xlsx(list(per_sample_editing = res),file = "first_look_edit_table.xlsx")
# 
# 
# # kd_test_tab <- kd_test_tab[all_variants < 4]
# kd_test_tab[,is_potential_edit := F]
# kd_test_tab[(ref == "A" & alt == "G") | (ref == "T" & alt == "C"),is_potential_edit := T]
# kd_test_tab[,is_AG := F]
# kd_test_tab[(ref == "A" & alt == "G") ,is_AG := T]
# kd_test_tab[,is_TC := F]
# kd_test_tab[(ref == "T" & alt == "C"),is_TC := T]
# 
# kd_test_tab[,.N,.(ref,alt)]
# 
# setorder(kd_test_tab,ref,alt)
# kd_test_tab[self_variant_count > 2][,.N,by = .(ref,alt,sample)]
# 
# edit_set <- kd_test_tab[(ref == "A" & alt == "G") | (ref == "T" & alt == "C")]
# control_set <- kd_test_tab[(ref == "G" & alt == "A") | (ref == "C" & alt == "T")]
# 
# control_set[,mean(self_variant_count),by = .(known_editting)]
# edit_set[,mean(self_variant_count),by = .(known_editting,)]
# 
# edit_set[self_variant_count > 2][,.N,by = sample]
# control_set[self_variant_count > 2][,.N,by = sample]
# 
# all_data_labels[,file := NULL]
# all_data_labels[,edit_count := sapply(res,nrow)]
# 
# res2 <- lapply(res,function(tab) {
#   tab <- tab[pop_AF < 0.02]
#   tab <- tab[AG_variants > 1 & all_variants < 4 | known_editting == T]
#   return(tab)
# })
# 
# 
# ###assign_cluster_var_size function
# # tab <- edit_set[sample == "u87_adar_ctrl__rep_1" & chr == "1" & ref == "A"]
# # pos <- tab$pos
# # 
# # assign_cluster_var_size <- function(pos,dist){
# #   sort_pos <- sort(pos)
# #   rle()
# #   diff(sort_pos)
# # }
# 
# 
# all_data_labels[,edit_count2 := sapply(res2,nrow)]
# 
# all_data_labels
# 
# table(tab$AG_variants)

