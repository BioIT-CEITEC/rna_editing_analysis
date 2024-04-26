library(data.table)
library(openxlsx)

#
# script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
# setwd(paste0(script_dir,"/../.."))


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
  fwrite(annot_tab,file = "editing_analysis/annot_tab.tsv")
  
  write.xlsx(per_sample_res,file = "editing_analysis/per_sample_editing.xlsx")
  write.xlsx(per_gene_res,file = "editing_analysis/per_coding_gene_editing.xlsx")
  write.xlsx(all_edit_sites_tab,file = "editing_analysis/all_editing_sites_tab.xlsx")

} 

know_editing_sites_analysis()


