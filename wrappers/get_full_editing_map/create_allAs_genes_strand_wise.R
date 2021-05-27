library(Biostrings)
library(GenomicRanges)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg38)
library(stringi)
genes_regions <- read.table("/mnt/ssd/ssd_3/references/homsap/GRCh38-p10/annot/GRCh38-p10_whole_genes.bed",stringsAsFactors = FALSE)

genes_regions <- genes_regions[genes_regions$V1 %in% c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y","MT"),]
genes_regions$V1 <- paste0("chr",genes_regions$V1)
genes_regions[genes_regions$V1 =="chrMT",]$V1 <- "chrM"
ranges <- GRanges(genes_regions$V1,IRanges(start=genes_regions$V2, end = genes_regions$V3),strand=genes_regions$V6)
seq_tab <- as.data.table(Views(Hsapiens,ranges))
seq_tab[seq_tab$strand == "-",]$dna <- reverse(seq_tab[seq_tab$strand == "-",]$dna)
as <- stringi::stri_locate_all(seq_tab$dna,regex = "[A]")
seq_tab$seqnames <- gsub("chr","",seq_tab$seqnames)
seq_tab$seqnames <- gsub("M","MT",seq_tab$seqnames)

seq_tab[,index := seq_along(seqnames)]

tab <- data.table(index = rep(seq_along(as),sapply(as,nrow)),rel_pos = unlist(lapply(as,function(x) x[,1])))
res <- merge(seq_tab[,.(seqnames,start,index)],tab,by = "index")
res[,gen_pos := start+rel_pos-1]
                  
final_pos_tab <- res[,c(2,5)]
write.table(final_pos_tab,"/mnt/ssd/ssd_3/references/homsap/GRCh38-p10/intervals/RNA/all_As_genes_strand_wise.pos",row.names = F,col.names = F,quote = F)

