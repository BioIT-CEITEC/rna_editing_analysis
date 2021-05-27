

load_full_editing_maps <- function(file_vec){
  all_per_pos <- list()
  all_per_read <- list()
  
  for(file in file_vec){
    sample_name <- gsub(".*/(.*).edit_tab.Rdata","\\1",file)
    load(file)
    per_pos_tab[,sample := sample_name]
    all_per_pos[[sample_name]] <- per_pos_tab
    
    per_read_tab[,sample := sample_name]
    all_per_read[[sample_name]] <- per_read_tab
  }
  
  all_per_pos <- rbindlist(all_per_pos)
  all_per_read <- rbindlist(all_per_read)
  
  return(list(per_pos = all_per_pos,per_read = all_per_read))
}
  
get_proximity_count <- function(vec,w = 500){
  tab <- data.table(c(vec,vec - w),c(seq_along(vec),seq_along(vec)))
  setorder(tab,V1)
  minus <- sapply(unique(tab$V2),function(x) sum(tab$V2[which(tab$V2 == x)[1]:which(tab$V2 == x)[2]] < x))
  tab <- data.table(c(vec,vec + w),c(seq_along(vec),seq_along(vec)))
  setorder(tab,V1)
  plus <- sapply(unique(tab$V2),function(x) sum(tab$V2[which(tab$V2 == x)[1]:which(tab$V2 == x)[2]] > x))
  res <- minus + plus
}

get_highly_modified_regions <- function(vec,window = 500,min_var_count = 5){
  
   counts <- get_proximity_count(vec,window)
   starts <- which(counts >= min_var_count & c(0,head(counts,-1)) < min_var_count)
   ends <- which(counts >= min_var_count & c(tail(counts,-1),0) < min_var_count)
   variant_count <- ends - starts +1 
   
   return(data.table(start_pos = vec[starts],end_pos = vec[ends],variant_count = variant_count))
}

run_all <- function(){
  file_vec <- c("/mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/first/input_files/results/u87_control_rep1.edit_tab.Rdata",
                "/mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/first/input_files/results/u87_control_rep2.edit_tab.Rdata",
                "/mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/first/input_files/results/u87_adar_kd_rep1.edit_tab.Rdata",
                "/mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/first/input_files/results/u87_adar_kd_rep2.edit_tab.Rdata")
  
  all_data <- load_full_editing_maps(file_vec)
  
  #
  #load Alu and snps and filter
  #
  SNP_tab <- fread("/mnt/ssd/ssd_3/references/homsap/GRCh38-p10/seq/GRCh38-p10.SNPs.bed.gz")
  all_data$per_pos <- merge(all_data$per_pos,SNP_tab[,.(chrom = V1,pos = V2,is_snp = TRUE)],by = c("chrom","pos"),all.x = T)
  all_data$per_pos[is.na(is_snp),is_snp := F]
  all_data$per_pos[is_snp == T,.N,by = sample]
  alu_regions <- fread("/mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/first/input_files/alu_pos.bed")
  alu_regions <- alu_regions[grep("^chr..?$",V1)]
  alu_regions[,V1 := gsub("^chr","",V1)]
  all_data$per_pos <- annotate_with_intervals(all_data$per_pos,alu_regions,"V4")
  setnames(all_data$per_pos,"V4","alu_region")
  all_data$per_pos[is.na(alu_region),alu_region := "no"]
  all_data$per_pos[alu_region != "no"]
  
  region_tab <- all_data$per_pos[,get_highly_modified_regions(pos,50,4),by = c("sample","chrom")]
  
  all_data$per_pos[,proximity_count_75 := get_proximity_count(pos,75),by = c("sample","chrom")]
  res_tab1 <- lapply(unique(all_data$per_pos$sample),function(x) table(all_data$per_pos[sample == x]$proximity_count_75))
  names(res_tab1) <- unique(all_data$per_pos$sample)
  res_tab1
  
  filtered_tab <- copy(all_data$per_read)
  filtered_tab <- merge(filtered_tab,all_data$per_pos[alu_region != "no"][,.(chrom,pos)],by = c("chrom","pos"))
  filtered_tab <- unique(filtered_tab,by = c("read","chrom","pos","sample"))
  filtered_tab[,cond := gsub("^.*?_(.*)_.*$","\\1",sample)]
  filtered_tab[,base_cap := toupper(base)]
  filtered_tab[,nuc_pop := .N,by = c("read","base_cap")]
  filtered_tab[,most_pop_nuc := base_cap[which.max(nuc_pop)],by = read]
  filtered_tab <- filtered_tab[base_cap == most_pop_nuc,]
  filtered_tab[,in_read_count := .N,by = c("sample","read")]

  
  test_filtering <- function(tab){
    print(nrow(tab))
    res <- tab[,.N,by = cond][,N_rel := N / sum(N)]
    print(res)
  }
  
  test_filtering(filtered_tab[in_read_count > 2 & in_read_count / read_mm > 0.5])
  test_filtering(filtered_tab)
  test_filtering(filtered_tab[in_read_count / read_mm > 0.2])
  test_filtering(filtered_tab[in_read_count > 1 &in_read_count / read_mm > 0.8])
  test_filtering(filtered_tab[in_read_count > 1])
  
  
  filtered_tab <- filtered_tab[in_read_count / read_mm > 0.2]
  res <- filtered_tab[,.N,by = cond][,N_rel := N / sum(N)]
  res
  res

  #split regions of intron overlaping reads
  
  regions <- filtered_tab[,.(chrom = chrom[1]
                            ,start = min(pos)
                            ,end = max(pos)
                            ,base = base_cap[1]
                            ,in_read_count = in_read_count[1]
                            ,read_mm = read_mm[1]
                            ,sample = sample[1]
                            ,cond = cond[1]),by = read]
  
  ranges <- GRanges( regions$chrom,IRanges(start=regions$start, end = regions$end),strand="+" )
  g <- as.data.table(reduce(ranges))
  
  res_tab2 <- lapply(unique(all_data$per_read$sample),function(x) table(all_data$per_read[sample == x]$in_read_count))
  names(res_tab2) <- unique(all_data$per_read$sample)
  res_tab2

  bam_vec <- c("/mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/first/input_files/mapped/u87_control_rep1.bam",
                "/mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/first/input_files/mapped/u87_control_rep2.bam",
                "/mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/first/input_files/mapped/u87_adar_kd_rep1.bam",
                "/mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/first/input_files/mapped/u87_adar_kd_rep2.bam")
  
  
  tab_orig <- as.data.table(read.xlsx("/mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/first/input_files/Supp_Tables_Erez.xlsx",sheet = 8))
  tab[,chr := trimws(chr)]
  tab
  
  tab <- fread("/mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/first/input_files/RNA_edit_erez_HG38.bed")
  tab <- tab[-156]
  pos_tab <- rbindlist(list(tab[1:196,.(chr = V1,pos = V2,pos_type = "r1_start",seq_id = 1:196)],
                            tab[1:196,.(chr = V1,pos = V3,pos_type = "r1_end",seq_id = 1:196)],
                            tab[197:392,.(chr = V1,pos = V2,pos_type = "r2_start",seq_id = 1:196)],
                            tab[197:392,.(chr = V1,pos = V3,pos_type = "r2_end",seq_id = 1:196)]))
  
  
  pos_tab[,chr := gsub("^chr","",chr)]
  gtf <- rtracklayer::import('/mnt/ssd/ssd_3/references/homsap/GRCh38-p10/annot/GRCh38-p10.gtf')
  gtf_df <- as.data.table(gtf)
  setnames(gtf_df,"type","gene_type")
  
  annot_tab <- annotate_with_intervals(pos_tab,gtf_df[gene_type == "gene"],"gene_type")
  annot_tab2 <- annotate_with_intervals(pos_tab,gtf_df[gene_type == "exon"],"gene_type")
  annot_tab2 <- annotate_with_intervals(annot_tab2,gtf_df[gene_type == "gene"],"gene_name")
  intron_tab <- annot_tab2[is.na(gene_type)][,gene_type := NULL]
  intron_tab <- dcast.data.table(intron_tab,formula = chr + seq_id + gene_name ~ pos_type,fill = 0,value.var = "pos")
  intron_tab <- intron_tab[r1_end != 0 & r2_end != 0 & r1_start != 0 & r2_start != 0]
  intron_tab <- rbind(intron_tab[,.(chr,start = r1_start,end = r1_end,gene_name,seq = "r1")],
                      intron_tab[,.(chr,start = r2_start,end = r2_end,gene_name,seq = "r2")])
  
  
  write.table(intron_tab,"/mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/first/input_files/intron_tab.bed",sep = "\t",row.names = F,col.names = F,quote = F)
  paste0("bedtools multicov -bams ",paste(bam_vec,collapse = " ")," -bed /mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/first/input_files/intron_tab.bed > /mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/first/input_files/coverage.tsv")
  
  gtf_to_bed <- gtf_df[gene_type == "gene" & seqnames %in% c(1:22,"X","Y"),1:12,with = F]
  write.table(gtf_to_bed,"/mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/first/input_files/all_gtf.bed",sep = "\t",row.names = F,col.names = F,quote = F)
  paste0("bedtools multicov -bams ",paste(bam_vec,collapse = " ")," -bed /mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/first/input_files/all_gtf.bed > /mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/first/input_files/gtf_coverage.tsv")
  gtf_cov_tab <- fread("/mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/first/input_files/gtf_coverage.tsv")
  # setnames(cov_tab,c("chr","start","end","seq","u87_control_rep1","u87_control_rep2","u87_adar_kd_rep1","u87_adar_kd_rep2"))
  gtf_mean_over_100 <- gtf_cov_tab[(V13 + V14 + V15 + V16) / 4 > 100]
  gtf_mean_over_100 <- gtf_cov_tab[V12 %in% unique(mean_over_50$seq)]
  t.test(c(gtf_mean_over_100$V13,gtf_mean_over_100$V14),c(gtf_mean_over_100$V15,gtf_mean_over_100$V16))
  
  
  cov_tab <- fread("/mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/first/input_files/coverage.tsv")
  setnames(cov_tab,c("chr","start","end","seq","gene_name","u87_control_rep1","u87_control_rep2","u87_adar_kd_rep1","u87_adar_kd_rep2"))
  mean_over_50 <- cov_tab[(u87_control_rep1 + u87_control_rep2 + u87_adar_kd_rep1 + u87_adar_kd_rep2) / 4 > 50]
  
  write.xlsx(list(all = cov_tab,mean_over_50 = cov_tab[(u87_control_rep1 + u87_control_rep2 + u87_adar_kd_rep1 + u87_adar_kd_rep2) / 4 > 50]),file = "/mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/first/input_files/coverage.xlsx")
  
  t.test(c(mean_over_50$u87_control_rep1,mean_over_50$u87_control_rep2),c(mean_over_50$u87_adar_kd_rep1,mean_over_50$u87_adar_kd_rep2))
  
  
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
  
  tab <- fread("/mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/first/input_files/mapped/u87_adar_kd_rep1.bam")
  tab[,V1 := gsub("chr","",V1)]
  
  tabs = apply(tab,1,function(x) all_data$per_pos[chrom == x[1] & pos > x[2] & pos < x[3],] )
  
  ranges <- GRanges( tab$chr,IRanges(start=tab$start1, end = tab$end1),strand="+" )
  seqs1 <- as.data.table(Views(Hsapiens,ranges))
  
  ranges <- GRanges( tab$chr,IRanges(start=tab$start2, end = tab$end2),strand="-" )
  seqs2 <- as.data.table(Views(Hsapiens,ranges))  
  
  
  library(data.table)
  library(openxlsx)
  
  gtf <- as.data.table(rtracklayer::import('/mnt/ssd/ssd_3/references/homsap/GRCh38-p10/annot/GRCh38-p10.gtf'))
  coverage <- as.data.table(read.xlsx("/mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/first/input_files/"))
  # gtf_gene <- gtf[type == "gene"]
  # gene_list <- apply(coverage,1, function(x) gtf_gene[x["chr"] == seqnames & as.numeric(x["start"]) > start &  as.numeric(x["start"]) < end]$gene_name)
  
 
  cov_filt <- coverage[(u87_control_rep1 + u87_control_rep2 + u87_adar_kd_rep1 + u87_adar_kd_rep2) / 4 > 50]
  setkey(cov_filt,chr,start,end)
  gene_bed <- fread("/mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/first/input_files/intron_tab.bed")
  setkey(gene_bed,V1,V2,V3)
  
  
  gene_bed <-  gene_bed[cov_filt]
  gene_bed <- gene_bed[1:8,]
  gene_bed[,exon_before := apply(gene_bed,1,function(x) max(as.numeric(gtf[x["V4"] == gene_name & as.numeric(x["V2"]) >  end  ]$exon_number),na.rm = T))]
  
  exon_bed <- gtf[type == "exon" & gene_name %in% unique(gene_bed$V4),.(seqnames,start,end,gene_name,exon_number,transcript_id)]
  exon_bed[,exon_length := .N + runif(1),by = transcript_id]
  exon_bed[,select := exon_length == max(exon_length),by = gene_name]
  exon_bed <- exon_bed[select == T]
  exon_bed[,select := NULL]
  exon_bed[,exon_length := NULL]
  
  write.table(exon_bed,"/mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/first/input_files/exon_bed.bed",sep = "\t",row.names = F,col.names = F,quote = F)
  paste0("bedtools multicov -bams ",paste(bam_vec,collapse = " ")," -bed /mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/first/input_files/exon_bed.bed > /mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/first/input_files/exon_bed_cov.tsv")
  exon_bed_cov_tab <- fread("/mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/RNA_edit_pipeline_test/first/input_files/exon_bed_cov.tsv")
  setnames(exon_bed_cov_tab,c("chr","start","end","gene_name","exon_num","trans","u87_control_rep1","u87_control_rep2","u87_adar_kd_rep1","u87_adar_kd_rep2")) 
  melt_exon_bed_cov_tab <- melt.data.table(exon_bed_cov_tab,measure.vars = c("u87_control_rep1","u87_control_rep2","u87_adar_kd_rep1","u87_adar_kd_rep2"),variable.name = "sample")
  
  melt_exon_bed_cov_tab[,norm_value := value / width * 100]
  ggplot(melt_exon_bed_cov_tab[gene_name == "CORO1C"], aes(x=exon_num, y=norm_value, group=sample)) +
    geom_line(aes(color=sample))+
    geom_point(aes(color=sample))
  
  
}









