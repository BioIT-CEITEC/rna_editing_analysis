library(data.table)


run_all <- function(args){
  out_filename <- args[1]
  panel_intervals_filename <- args[2]
  panel_snps_filename <- args[3]
  panel_cnv_data_filename <- args[4]
  cov_tab_filenames <- args[5:(which(args == "snps") - 1)]
  snp_tab_filenames <- args[(which(args == "snps") + 1):length(args)]
  
  #additional args
  trisomy_prob=0.01
  quad_prob=0.005
  
  panel_intervals <- fread(panel_intervals_filename)
  setnames(panel_intervals,c("V1","V2","V3","V4"),c("chr","start","end","region_name"))
  panel_intervals[,V5 := NULL]
  panel_intervals[,V6 := NULL]
  panel_intervals[,region_name_suffix := paste0("_",as.character(seq_along(start))),by = region_name]
  panel_intervals[,region_name_count := .N,by = region_name]
  panel_intervals[region_name_count == 1,region_name_suffix := ""]
  panel_intervals[,region_name := paste0(region_name,region_name_suffix)]
  panel_intervals[,c("region_name_count","region_name_suffix") := NULL]
  setkeyv(panel_intervals, c("chr","start","end"))
  
  snp_tab <- fread_vector_of_files(snp_tab_filenames,".*\\/(.*)_snpAF.tsv")
  setnames(snp_tab,c("chrom","pos","A","C","G","T","cov","sample"))
  snp_tab <- snp_tab[cov > 0]
  panel_snps <- fread(panel_snps_filename)
  setnames(panel_snps,c("chrom","pos","ref","alt","AF"))
  snp_tab <- merge(snp_tab,panel_snps,by = c("chrom","pos"))
  snp_tab[,alt_count := as.numeric(apply(snp_tab,1,function(x) x[x["alt"]]))]
  snp_tab[,VAF := alt_count / cov]
  
  #genetics law
  snp_tab[,pop_HET_probability := 2*AF*(1-AF)]
  #AF = 0,25,33,50,66,75,100
  prob_vec <- c(0,25,33,50,66,75,100)/100
  prior_prob <- matrix(0,ncol = 7,nrow = nrow(snp_tab))
  prior_prob[,1] <- (1 - snp_tab$pop_HET_probability) / 2
  prior_prob[,7] <- prior_prob[,1]
  prior_prob[,2] <- (snp_tab$pop_HET_probability * quad_prob) / 3
  prior_prob[,6] <- prior_prob[,2]
  prior_prob[,3] <- (snp_tab$pop_HET_probability * trisomy_prob) / 2
  prior_prob[,5] <- prior_prob[,3]
  prior_prob[,4] <- snp_tab$pop_HET_probability - 2*prior_prob[,2] - 2*prior_prob[,3]
  
  pdata <- sapply(prob_vec,function(prob) {
      mapply(function(alt_count,cov) prob^alt_count * (1 - prob)^(cov - alt_count),alt_count = snp_tab$alt_count,cov = snp_tab$cov)
  })
  
  p_result <- prior_prob * pdata
  p_result <- p_result / rowSums(p_result)
  snp_tab[,HOM_prob := p_result[,1] + p_result[,7]]
  snp_tab[,HET_prob := p_result[,4] - (prior_prob[,2] + prior_prob[,6])]
  snp_tab[,TRI_prob := p_result[,3] + p_result[,5]]
  snp_tab[,QUAD_prob := 2*(prior_prob[,2] + prior_prob[,6])]

  snp_tab[,pos2 := pos]
  snp_tab <- foverlaps(snp_tab, panel_intervals, by.x=c("chrom","pos","pos2"), by.y=c("chr","start","end"), nomatch = 0)
  snp_tab[,pos2 := NULL]
  
  print("snp_end")
  timestamp()
  
  #? GC content fit
  # genome=BSgenome.Hsapiens.UCSC.hg19
  # panel_ranges <- GRanges( paste0("chr",panel_intervals$chr),
  #                           IRanges(start=panel_intervals$start, end = panel_intervals$end),
  #                           strand="+" )
  # 
  # af <- alphabetFrequency(Views(genome,panel_ranges), baseOnly = TRUE, as.prob = TRUE)
  # panel_intervals[,GC_content := af[,"C"] + af[,"G"]]
  # rm(af)
  
  cov_tab <- fread_vector_of_files(cov_tab_filenames,".*\\/(.*)_cov.tsv")
  setnames(cov_tab,c("V1","V2","V3"),c("chr","pos","cov"))
  cov_tab[,pos2 := pos]
  cov_tab <- foverlaps(cov_tab, panel_intervals, by.x=c("chr","pos","pos2"), by.y=c("chr","start","end"), nomatch = 0)
  # cov_tab[,middle_pos := pos[round(.N/2)],by = .(chr,start,end,sample,region_name)]
  # estimate_region_coverage <- function(region_cov_vec,region_pos_vec){
  #   # https://stats.stackexchange.com/questions/83022/how-to-fit-data-that-looks-like-a-gaussian
  #   FIT <- fitdist(rep(region_pos_vec - min(region_pos_vec),times = region_cov_vec - min(region_cov_vec)), "norm")
  #   cov <- mean(region_cov_vec[which(round(FIT$estimate[1]) == region_pos_vec) + -2:2])
  #   return(cov)
  # }
  region_cov_tab <- cov_tab[,.(cov = ifelse(.N > 10,mean(cov[round(.N/2) + -4:4]),mean(cov))),by = .(sample,region_name)]
  setkey(region_cov_tab)
  rm(cov_tab)
  region_cov_tab <- region_cov_tab[CJ(unique(region_cov_tab$sample),unique(region_cov_tab$region_name))] 
  region_cov_tab[is.na(cov), cov := 0]
  region_cov_tab[,cov_raw := cov]
  overall_mean <- mean(region_cov_tab$cov_raw)
  region_cov_tab[,cov := cov_raw / mean(cov_raw) * overall_mean,by = sample]
  
  
  #per number of clusters
  # for testing lets assumme k=3 in the future test values 1-4 ?how to test?
  # k_clusters = 1
  
  # region_meancov_fit <- fitdist(round(region_cov_tab$cov),"nbinom")
  # credible_values_threshold <- qnbinom(c(0.05,0.95),mu = region_meancov_fit$estimate["mu"],size = region_meancov_fit$estimate["size"])
  # region_cov_tab[,is_credible_region := all(cov > credible_values_threshold[1] & cov < credible_values_threshold[2]),by = region_name]
  
  # #fuzzy clustering of PCA
  # mat <- as.matrix(dcast.data.table(region_cov_tab[is_credible_region == T,],formula = region_name ~ sample,value.var = "cov",fill = 0)[,-1,with = F])
  # res <- prcomp(t(mat),scale. = T)
  # clustered_object <- cluster::fanny(res$x,k_clusters)
  # autoplot(clustered_object,frame = T)
  
  # hierarchical clustering on correlation distances
  mat <- as.matrix(dcast.data.table(region_cov_tab,formula = region_name ~ sample,value.var = "cov",fill = 0)[,-1,with = F])
  hclust <- hclust(-simil(mat,method = cor,by_rows = F,use = "complete.obs"))
  clustering <- cutree(hclust,h = -0.90)
  k_clusters <- length(unique(clustering))
  
  # GC content does not corelate
  # setkey(panel_intervals,region_name)
  # panel_intervals[,est_cov := GC_content / sum(GC_content) * overall_mean]
  # CG_clust <- simil(t(mat),t(panel_intervals[unique(region_cov_tab$region_name)]$est_cov),method = cor)[,1]
  
  # sample_cluster_membership <- as.data.table(clustered_object$membership,keep.rownames = T)
  # setnames(sample_cluster_membership,c("sample",paste("clust",1:k_clusters,sep = "_")))
  # setkey(sample_cluster_membership,sample)
  # # sample_cluster_membership[,clustering := paste("clust",clustered_object$clustering[order(names(clustered_object$clustering))],sep = "_")]
  # setkey(region_cov_tab,region_name,sample)
  # setkey(region_cov_tab,region_name)
  # test_mat <- mat[10:17,1:10]
  # 
  # pheatmap(round(test_mat, display_numbers = T, color = colorRampPalette(c('white','red'))(100), cluster_rows = F, cluster_cols = F, fontsize_number = 15)
  
  
  if(k_clusters > 1){
    distance_faktor <- 8
    sim_matrix <- as.matrix(simil(mat,method = cor,by_rows = F,diag = T,upper = T,use = "complete.obs"))
    diag(sim_matrix) <- 1
    sim_matrix <- apply(sim_matrix,1,function(x) (x - 0.5) * 2)
    sim_matrix[sim_matrix < 0] <- 0
    sim_matrix <- apply(sim_matrix,1,function(x) x ^ distance_faktor)
    sim_matrix <- sim_matrix / rowSums(sim_matrix)
    sim_matrix <- round(sim_matrix * 1000)
    
    get_one_region_sample_clust_nbinom_fit_coefs <- function(coverage_vec,weight_vec){
      return(as.list(fitdist(rep(round(coverage_vec),weight_vec),"nbinom")$estimate))
    }
    
    get_fit_dist_estimates <- function(sample_vec,cov_vec,sim_matrix,by = NULL){
      print(by)
      fit_dist_estimates <- mclapply(sample_vec,function(x) get_one_region_sample_clust_nbinom_fit_coefs(cov_vec,sim_matrix[x,]),
               mc.preschedule = TRUE, mc.set.seed = TRUE,
               mc.silent = TRUE, mc.cores = 20,
               mc.cleanup = TRUE)
      fit_dist_estimates <- simplify2array(fit_dist_estimates)
      if(length(dim(fit_dist_estimates)) == 2 && dim(fit_dist_estimates)[1] == 2){
        return(list(size = fit_dist_estimates[1,],mu = fit_dist_estimates[2,]))
      } else {
        return(list(size = NA,mu = NA))
      }
    }

    print("start_fitting")
    timestamp()
    region_cov_tab[,c("size","mu") := get_fit_dist_estimates(sample,cov,sim_matrix,.BY[[1]]),by = region_name]
    copies_prob_tab <- region_cov_tab
    print("end_fitting")
    timestamp()
    # get_one_region_sample_clust_nbinom_fit_coefs <- function(coverage_vec,weight_vec){
    #   return(as.list(fitdist(rep(round(coverage_vec),round(weight_vec * 100)),"nbinom")$estimate))
    # }
    # 
    # get_fit_dist_estimates <- function(cluster_id,sample_cluster_membership,region_cov_tab){
    #   weight_vec <- sample_cluster_membership[[cluster_id]]
    #   fit_dist_estimates <- region_cov_tab[,get_one_region_sample_clust_nbinom_fit_coefs(cov,weight_vec),by = region_name]
    #   return(fit_dist_estimates)
    # }
    # 
    # all_clust_fit_dist_estimates <- lapply(names(sample_cluster_membership)[-1],get_fit_dist_estimates,sample_cluster_membership,region_cov_tab)
    # names(all_clust_fit_dist_estimates) <- names(sample_cluster_membership)[-1]
    # all_clust_fit_dist_estimates <- rbindlist(all_clust_fit_dist_estimates,use.names = T,idcol = "cluster")
    # setkey(all_clust_fit_dist_estimates,region_name)
    # 
    # region_cov_tab <- region_cov_tab[all_clust_fit_dist_estimates,allow.cartesian=TRUE]
    # sample_cluster_membership <- melt.data.table(sample_cluster_membership,id.vars = "sample",variable.name = "cluster",value.name = "membership")
    # region_cov_tab <- merge(region_cov_tab,sample_cluster_membership,by = c("sample","cluster"))
    # copies_prob_tab <- region_cov_tab[,.(cov = cov[1],mu = sum(mu * membership),size = sum(size * membership)),by = .(sample,region_name)]
  } else {
    region_cov_tab[,c("size","mu") := as.list(fitdist(round(cov),"nbinom")$estimate),by = region_name]
    copies_prob_tab <- region_cov_tab[,c("is_credible_region","cov_raw") := NULL]
  }

  copies_prob_tab[,copies_0 := dnbinom(x = round(cov),mu = 1,size = 0.1)]
  copies_prob_tab[,copies_1 := dnbinom(x = round(cov),mu = mu / 2 * 1,size = size / (2 / 1)^2 )]
  copies_prob_tab[,copies_2 := dnbinom(x = round(cov),mu = mu / 2 * 2,size = size / (2 / 2)^2 )]
  copies_prob_tab[,copies_3 := dnbinom(x = round(cov),mu = mu / 2 * 3,size = size / (2 / 3)^2 )]
  copies_prob_tab[,copies_4 := dnbinom(x = round(cov),mu = mu / 2 * 4,size = size / (2 / 4)^2 )]
  
  
  copies_prob_tab <- cbind(copies_prob_tab[,grep("copies",names(copies_prob_tab),invert = T,value = T),with = F],
                           t(apply(copies_prob_tab[,grep("copies",names(copies_prob_tab),value = T),with = F],1,function(x) x / sum(x))))
  
  
  save(copies_prob_tab,snp_tab,panel_intervals,file = "stage277_BRONCO_CNV.test1/CNV_panel_analysis/results/test_res.R")
  # copies_prob_tab
  # 
  # 
  # x <- data.table(copies_1=rnbinom(10000,mu = mu / 2 * 1,size = size / (2 / 1)^2 ),
  #                 copies_2=rnbinom(10000,mu = mu / 2 * 2,size = size / (2 / 2)^2 ),
  #                 copies_3=rnbinom(10000,mu = mu / 2 * 3,size = size / (2 / 3)^2 ),
  #                 copies_4=rnbinom(10000,mu = mu / 2 * 4,size = size / (2 / 4)^2 ))
  # data <- melt(x)
  # ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.25)
  # 
  # # one_sample_clust_nbinom_fit_coefs <- function(cluster_id,region_cov_tab,sample_cluster_membership){
  # #   weight_vec <- sample_cluster_membership[[cluster_id]]
  # #   fit_dist_estimates <- region_cov_tab[,get_one_region_sample_clust_nbinom_fit_coefs(cov,weight_vec),by = region_name]
  # #   return(fit_dist_estimates)
  # # }
  # 
  # # one_sample_clust_nbinom_fit_coefs <- function(cluster_id,region_cov_tab,sample_cluster_membership){
  # #   fit_dist_estimates <- region_cov_tab[sample %in% sample_cluster_membership[clustering == cluster_id]$sample,as.list(fitdist(round(cov),"nbinom")$estimate),by = region_name]
  # #   return(fit_dist_estimates)
  # # }
  # 
  # per_clust_fit_coefs <- lapply(sort(unique(sample_cluster_membership$clustering)),one_sample_clust_nbinom_fit_coefs,region_cov_tab,sample_cluster_membership)
  # names(per_clust_fit_coefs) <- sort(unique(sample_cluster_membership$clustering))
  # per_clust_fit_coefs <- rbindlist(per_clust_fit_coefs,use.names = T,idcol = "clustering")
  # 
  # 
  # 
  # for(cluster_number in 1:k_clusters){
  #   fit_dist_estimates <- region_cov_tab[,get_one_region_sample_clust_nbinom_fit_coefs(cov,sample_cluster_membership[sample,][[paste0("V",cluster_number)]]),by = region_name]
  #   region_cov_tab[sample %in% names(clustered_object$clustering)[clustered_object$clustering == cluster_number]]
  #   
  # }
  # 
  # # 
  # # region_cov_vec <- cov_tab[region_name == "BRCC3_21" & sample == "BR-0737"]$cov
  # # region_pos_vec <- cov_tab[region_name == "BRCC3_21" & sample == "BR-0737"]$pos
  # 
  # 
  # # estimate_region_coverage <- function(region_cov_vec,region_pos_vec){
  # #   # https://stats.stackexchange.com/questions/83022/how-to-fit-data-that-looks-like-a-gaussian
  # #   FIT <- fitdist(rep(region_pos_vec - min(region_pos_vec),times = region_cov_vec - min(region_cov_vec)), "norm")
  # #   cov <- mean(region_cov_vec[which(round(FIT$estimate[1]) == region_pos_vec) + -2:2])
  # #   return(cov)
  # # }
  # 
  # 
  # cnv_intervals <- rbindlist(lapply(unique(cov_tab$sample),function(x){
  #   tab <- merge(cnv_intervals,snp_tab[sample == x,.(chr = V1,snp_pos = V2,VAF,snp_cov = (ref + alt))],by = c("chr","snp_pos"),all.x = T)
  #   tab <- merge(tab,cov_tab[sample == x,.(reg_name,mean_cov,cov_sd,max_cov)],by = c("reg_name"),all.x = T)
  #   tab[,sample := x]
  #   return(tab)
  # }))
  # 
  # cnv_intervals[is.na(VAF),c("VAF","snp_cov") := 0]
  # cnv_intervals[is.na(mean_cov),c("mean_cov","cov_sd","max_cov") := 0]
  # 
  # sample_size_norm_tab <- cnv_intervals[,.(norm_factor = sum(mean_cov)),by =sample]
  # sample_size_norm_tab[,norm_factor := norm_factor / mean(norm_factor)]
  # cnv_intervals <- merge(cnv_intervals,sample_size_norm_tab,by = "sample")
  # cnv_intervals[,norm_mean_cov := mean_cov / norm_factor]
  # cnv_intervals[,norm_max_cov := max_cov / norm_factor]
  # # cnv_intervals[,norm_factor := NULL]
  # 
  # cnv_intervals <- merge(cnv_intervals,cnv_intervals[chr %in% c("Y"),.(is_woman = median(mean_cov) < 1),by = sample],by = "sample")
  # 
  # cnv_intervals[!(chr %in% c("X","Y")),copy_estimate := norm_max_cov / mean(norm_max_cov) * 2,by = reg_name]
  # cnv_intervals[chr == "X",copy_estimate := norm_max_cov / mean(norm_max_cov) * (as.integer(is_woman) + 1),by = c("reg_name","is_woman")]
  # cnv_intervals[chr == "Y",copy_estimate := norm_max_cov / mean(norm_max_cov),by = c("reg_name")]
  # cnv_intervals[,copy_estimate_rollmean_50 := rollmean(copy_estimate,k = 50,fill = 0),by = c("sample","chr")]
  # 
  # cnv_intervals[,het_snp := F]
  # cnv_intervals[snp_cov > 5,het_snp := VAF < 0.75 & VAF > 0.25]
  # cnv_intervals[,het_sum_50 := rollsum(as.integer(het_snp),k = 50,fill = 0),by = c("sample","chr")]
  # 
  # dir.create(dirname(out_filename),recursive = T)
  # save(cnv_intervals,file = out_filename)

}


# develop and test
args <- c("stage277_BRONCO_CNV.test1/CNV_panel_analysis/results/CNV_result_tab.tsv","/mnt/ssd/ssd_3/references/homsap/GRCh37-p13/intervals/BRONCO/BRONCO.bed","/mnt/ssd/ssd_3/references/homsap/GRCh37-p13/other/snp/BRONCO/BRONCO_snps.bed","/mnt/ssd/ssd_3/references/homsap/GRCh37-p13/other/cnv_panel_data/BRONCO/BRONCO_cnv_panel_data.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/H-0251_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/10239_19_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0686_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0762_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/H-0224_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/H-0233_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/H-0225_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/H-0244_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/H-0227_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0760_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0778_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/H-0233-BS_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0708_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/H-0229_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0733_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0695_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/H-0230_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0723_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/H-0247_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0781_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0732_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/H-0239_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0696_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0699_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/H-0213_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0774_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0777_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/H-0231_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/PZ-3_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0689_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0737_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0712_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0782_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0721_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/H-0228_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0773_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0775_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0731_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0738_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0784_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/H-0222_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0714_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0709_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0783_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/H-0215_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0758_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/H-0232_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/H-0214_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0669_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0766_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0720_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0676_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0671_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0768_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0715_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0678_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0705_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/H-0243_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0674_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BB17709_20_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0698_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0772_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0767_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0776_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0736_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/H-0217_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0727_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0679_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0716_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/H-0219_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0670_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0729_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/971_19_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/H-0238_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0688_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0749_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/H-0223_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0734_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0703_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0687_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0763_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0697_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0713_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0787_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/H-0245_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/H-0246_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0735_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0726_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0742_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0764_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0769_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0785_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0725_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0748_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0759_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0672_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0673_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0710_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/ZB-2_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0719_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0675_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0765_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0706_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0707_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/MI-1_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0780_cov.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/depth_tab/BR-0779_cov.tsv","snps","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0764_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/H-0232_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/MI-1_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/H-0229_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/H-0225_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0766_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/H-0245_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0736_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0673_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BB17709_20_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0738_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/H-0238_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/H-0246_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/H-0230_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0767_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/H-0215_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0783_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0749_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0732_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/H-0247_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/H-0233_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0695_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0699_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0782_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0712_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0765_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/PZ-3_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0769_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0775_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/H-0228_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0729_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0669_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0727_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/H-0224_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/H-0214_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/H-0243_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0773_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0737_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0742_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0705_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0758_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/H-0251_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0776_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0703_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0707_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0720_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0689_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0675_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0760_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/971_19_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0735_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0697_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0731_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0779_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0748_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/H-0223_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0678_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0781_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0762_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/H-0239_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0713_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0787_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0763_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/ZB-2_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/H-0219_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0734_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0706_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0686_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0708_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0772_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/H-0227_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0723_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0777_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0688_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/H-0231_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0778_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0774_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0709_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0733_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0671_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/H-0213_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0784_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0679_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0719_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/H-0217_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0670_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/H-0233-BS_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0714_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0725_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0710_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0726_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0715_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0716_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0676_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0721_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0785_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/H-0244_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0687_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0696_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/10239_19_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0759_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0674_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0780_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0698_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0768_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/H-0222_snpAF.tsv","stage277_BRONCO_CNV.test1/CNV_panel_analysis/snp_tab/BR-0672_snpAF.tsv")
setwd("/mnt/ssd/ssd_1/snakemake/")
# args <- "/mnt/ssd/ssd_1/snakemake/stage172_test_CNV_analysis/CNV_analysis/results/CNV_tabs.Rdata"
# args <- c(args,"/mnt/ssd/ssd_3/references/homsap/GRCh38-p10/other/cnv_intervals/CNV_LOH.tsv")
# args <- c(args,list.files("/mnt/ssd/ssd_1/snakemake/stage172_test_CNV_analysis/CNV_analysis/depth_tab",pattern = ".tsv",full.names = T))
# args <- c(args,"snps")
# args <- c(args,list.files("/mnt/ssd/ssd_1/snakemake/stage172_test_CNV_analysis/CNV_analysis/snp_tab",pattern = ".tsv",full.names = T))
script_dir <- "/mnt/nfs/shared/999993-Bioda/bioda_snakemake/wraps/CNV_analysis/cnv_computation"

#run as Rscript

# script_dir <- dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))
# args <- commandArgs(trailingOnly = T)
# run_all(args)
# print("start")
# timestamp()
# run_all(args)
# print("end")
# timestamp()

# system("which R")
