#2022-05-05
#Plot amount of accessible chromatin and amount of genetic differences in ACRs

library(ggplot2)
library(RColorBrewer)
library(gplots)
library(tidyr)
library(ggforce)
library(ggpubr)
library(dplyr)
library(ggrepel)


ATAC_peak_files <- list.files("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-05-05_Downsample_ACR_analysis/DownsampledBams/Peaks", pattern = "peaks.txt",full.names = TRUE)

ATAC_peak_files <- ATAC_peak_files[-c(2,8)]



ALL <- lapply(ATAC_peak_files, function(i){
 read.delim(i, header =FALSE)
})


#adding read numbers from atac bam files manually
read_number <-  read.table("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Downsampled_acrs/2023-10-17_downsampled_bam_readcounts.txt", quote="\"", comment.char="")




samples <- read.table("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Downsampled_acrs/2023-10-17_downsampled_bam_filenames.txt", quote="\"", comment.char="")



Summary <- list()

for (i in 1:19){
 ALL[[i]]$Size <- ALL[[i]]$V3 - ALL[[i]]$V2
 ALL[[i]]$Genotype <- gsub("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-05-05_Downsample_ACR_analysis/DownsampledBams/Peaks/(.*?)", "\\1",ATAC_peak_files[[i]])
 ALL[[i]]$Genotype <- gsub(".05.peaks.txt","",ALL[[i]]$Genotype)
 ALL[[i]]$Number_of_Reads <-read_number[i,]
 ALL[[i]]$Open_chromatin_bp <- sum(ALL[[i]]$Size)
 ALL[[i]]$Number_of_Peaks <- dim(ALL[[i]])[1]
 Summary[[i]] <- unique(ALL[[i]][,c(12:15)])
 
 
}

All_ATAC_sum <- do.call(rbind, Summary)
All_ATAC_sum$Market_Type <- c("Mandarin","SourOrange","Papeda","Papeda","Lemon","Mandarin","Pummelo","Mandarin","Microcitrus","Mandarin","Grapefruit","Mandarin","Citron","Mandarin","Trifoliate","Mandarin","Mandarin","Grapefruit","Orange") 

All_ATAC_sum$Type <- factor(All_ATAC_sum$Market_Type, levels=c("Mandarin","Papeda","Citron","Grapefruit","Lemon","Pummelo","Trifoliate","Microcitrus","Orange"))

All_ATAC_sum$Method <- "Genrich_0.05"

pdf("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Summarize_ATAC/bp_accessible_chromatin.pdf")
ggplot(All_ATAC_sum, aes(x=reorder(Genotype,-Open_chromatin_bp), y = (Open_chromatin_bp/1000000), fill =Market_Type)) + geom_bar(stat="identity") + coord_flip() + theme_classic() + ylab("Accessible chromatin (Mb)") + xlab("Genotype") + theme(text=element_text(size=20)) + scale_fill_brewer(palette="Set3")

dev.off()



c <- ggplot(All_ATAC_sum, aes(x=Number_of_Reads/1000000, y = Open_chromatin_bp/1000000)) + geom_point(size=1)+ geom_label_repel(aes(label=Genotype)) + theme_bw() + ylab("bp open chromatin (Mb)") + xlab("Number of reads (Million reads)") + theme(text=element_text(size=20))

d <- ggplot(All_ATAC_sum, aes(x=Number_of_Reads/1000000, y = Open_chromatin_bp/1000000)) + geom_point(size=1) + theme_bw() + ylab("bp open chromatin (Mb)") + xlab("Number of reads (Million reads)") + theme(text=element_text(size=20))

pdf("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Summarize_ATAC/Corr_num.reads.openchrom.pdf")
c + stat_cor(method="spearman")

d + stat_cor(method="spearman")

dev.off()







