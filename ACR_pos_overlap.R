#!/usr/bin/env Rscript

#args[1] is strict ACR calling
#args[2] is loose ACR calling
#args[3] is MACS2 ACR calling

#parse args
args = commandArgs(trailingOnly=TRUE)


library(dplyr)
library(scales)
#library(karyoploteR)
library(ggplot2)

print("ACR OVERLAP")

genrich_loose <- read.table("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_by_geno/atac.occupancy.dist.txt")

pdf("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACR_OVERLAP/2023-09-05_ACR.pos.overlap.pdf")

plot(genrich_loose[,1], genrich_loose[,2] / sum(genrich_loose[,2]), 'h',col="darkred",lwd=4, xlab="No. of assayed genotypes",ylab="Fraction of bases",main="Socrate_Peaks")

dev.off()

Merged_ACRs.counts <- read.delim("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-05-05_Downsample_ACR_analysis/ACR_overlap/Merged_ACRs.counts.bed", header=FALSE)

Cat_all_acrs <- read.delim("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-05-05_Downsample_ACR_analysis/ACR_overlap/Cat_all_acrs.bed", header=FALSE)



Merged_ACRs.counts$Size <- Merged_ACRs.counts$V3 - Merged_ACRs.counts$V2


pdf("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-05-05_Downsample_ACR_analysis/ACR_overlap/2023-05-29_ACR.size.pdf")
#how does size of acrs change after collapsing overlapping ACRs.
Cat_all_acrs$Size <- Cat_all_acrs$V3 - Cat_all_acrs$V2
size_summary <- summary(Cat_all_acrs$Size)[4]
ggplot(Cat_all_acrs, aes(x=Size)) + geom_density() + theme_classic() + xlab("Size of all MERGED ACRs") + ggtitle(size_summary)

size_summary_merge <- summary(Merged_ACRs.counts$Size)[4]
ggplot(Merged_ACRs.counts, aes(x=Size)) + geom_density()  + theme_classic() + xlab("Size of all unmerged ACRs") + ggtitle(size_summary_merge)
summary(Merged_ACRs.counts$Size)

dev.off()


#make column denoting whether an collapsed ACR overlaps atleast one ACR in a specific genotype, 0 means no , 1 means yes
Merged_ACRs.counts$Overlap <- ifelse(Merged_ACRs.counts$V5 >= 1, 1, 0)



Overlap_summary <- Merged_ACRs.counts %>% group_by(V1,V2,V3) %>% summarize(Freq = sum(Overlap))

#This produces table of unique ACRs and the genotypes they belong to
solos <- Overlap_summary %>% filter(Freq == 1)
solos_counts <- merge(Merged_ACRs.counts,solos , by=c("V1","V2","V3")) %>% filter(V5==1)

solos_table <-as.data.frame(table(solos_counts$V4))

z <- solos_counts %>% filter(V4 == "Interdonato_lemon")
Others <- solos_counts %>% filter(V4 != "Interdonato_lemon")

cclem_genome <- read.delim("/bigdata/seymourlab/idiaz026/REFERENCES/Cclementina/Cclementine_genome.bed",header=F)[c(1:9),]

#custom.genome <- toGRanges(cclem_genome)
#kp <- plotKaryotype(genome = custom.genome)
#AS_ACRs <- toGRanges(z)
#all_acrs <- toGRanges(Others)

#kp <- plotKaryotype(genome=custom.genome, plot.type=2)

#kpPlotRegions(kp, data=AS_ACRs,col="red",data.panel=1)
#kpPlotRegions(kp, data=all_acrs,data.panel = 2)


pdf("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-05-05_Downsample_ACR_analysis/2023-05-29_singleton_ACRs_bygeno.pdf")
ggplot(solos_table, aes(x=reorder(Var1, -Freq ), y=Freq )) + geom_bar(stat="identity") + coord_flip()+ theme_bw() + xlab("Genotype") + ylab("# Unique ACRs")

dev.off()


Unique_ACRs <- Overlap_summary %>% filter(Freq ==1)
Common_ACRs <- Overlap_summary %>% filter(Freq > 1 & Freq < 17)
High_freq_ACRs <- Overlap_summary %>% filter(Freq >= 17)

Unique_ACRs$Type <- "Unique"
Common_ACRs$Type <- "Common"
High_freq_ACRs$Type <-  "High_frequency"

ACRs_classified <- rbind(Unique_ACRs,Common_ACRs,High_freq_ACRs)[,c(-4)]

write.table(ACRs_classified, "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-05-05_Downsample_ACR_analysis/ACRs_classified/ACRs_classified.bed" , sep ="\t", col.names=FALSE, row.names=FALSE,quote=FALSE)

write.table(Unique_ACRs[,c(1:3)], "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-05-05_Downsample_ACR_analysis/ACRs_classified/Unique_ACRs.bed", sep ="\t", col.names=FALSE, row.names=FALSE,quote=FALSE)

write.table(Common_ACRs[,c(1:3)], "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-05-05_Downsample_ACR_analysis/ACRs_classified/Common_ACRs.bed",sep="\t",col.names = FALSE, row.names=FALSE,quote =FALSE)

write.table(High_freq_ACRs[,c(1:3)], "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-05-05_Downsample_ACR_analysis/ACRs_classified/High_frequency_ACRs.bed",sep="\t",col.names = FALSE, row.names=FALSE,quote =FALSE)

png("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-05-05_Downsample_ACR_analysis/ACR_overlap/2023-05-29_ACR_frequency.png")
ggplot(Overlap_summary,aes(x=Freq)) + geom_histogram(bins = 39)+ theme_classic() + scale_x_continuous(breaks = scales::pretty_breaks(n = 17)) + xlab("ACR frequency across 20 genotypes") + theme(text=element_text(size=15))

dev.off()

head(Overlap_summary)

#ACR present in one individual = unique, 2 - 16 = common, > 16 = high freq
