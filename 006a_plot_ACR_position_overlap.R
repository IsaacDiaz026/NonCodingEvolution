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
library(UpSetR)
library(stringr)
library(GenomicRanges)

mutations <- read.csv( system.file("extdata", "mutations.csv", package = "UpSetR"), header=T, sep = ",")


genrich_loose <- read.table("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACR_OVERLAP/atac.occupancy.dist.txt")

pdf("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACR_OVERLAP/2023-09-07_ACR.pos.overlap.pdf")

plot(genrich_loose[,1], genrich_loose[,2] / sum(genrich_loose[,2]), 'h',col="darkred",lwd=4, xlab="No. of assayed genotypes",ylab="Fraction of bases",main="Genrich_0.05")

dev.off()

Merged_ACRs.counts <- read.delim("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACR_OVERLAP/Merged_ACRs.counts.bed", header=FALSE)

Cat_all_acrs <- read.delim("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACR_OVERLAP/Cat_all_acrs.bed", header=FALSE)


Merged_ACRs.counts$Size <- Merged_ACRs.counts$V3 - Merged_ACRs.counts$V2


#Filter ACRs that are clustered - this leads to ACR boundaries expanding which biases ACR overlap

b <- Merged_ACRs.counts %>% filter(V5 > 0)

hist(b$V5, xlab = "# of ACRs in one sample overlapping consensus ACR")

b_to_plot <- unique(b[,c(1,2,3,6)],)

c <- Merged_ACRs.counts %>% filter(V5 < 2)

c_to_plot <- unique(c[,c(1,2,3,6)],)

pdf("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACR_OVERLAP/2023-09-06_ACR.size.pdf")
#how does size of acrs change after collapsing overlapping ACRs.
Cat_all_acrs$Size <- Cat_all_acrs$V3 - Cat_all_acrs$V2
size_summary <- summary(Cat_all_acrs$Size)[3]
ggplot(Cat_all_acrs, aes(x=Size)) + geom_density() + theme_classic() + xlab("Size of all unmerged ACRs") + ggtitle(paste0("Median_ACR_size_",size_summary," n=",nrow(Cat_all_acrs), sep=""))+ xlim(0,2000)


size_summary_merge <- summary(b_to_plot$Size)[3]
ggplot(b_to_plot, aes(x=Size)) + geom_density()  + theme_classic() + xlab("Size of all consensus ACRs") + ggtitle(paste0("Median_ACR_size_",size_summary_merge," n=",nrow(b_to_plot), sep=""))+ xlim(0,2000)
summary(b_to_plot$Size)


num_filt <- dim(b_to_plot)[1] - dim(c_to_plot)[1]
hist(b$V5, xlab = "# of ACRs in one sample overlapping consensus ACR",main =paste0(num_filt," consensus ACRs filtered for being clustered",sep="" ))

size_summary_filt <- summary(c$Size)[3]
ggplot(c_to_plot, aes(x=Size)) + geom_density()  + theme_classic() + xlab("Size of all consensus, filtered ACRs") + ggtitle(paste0("Median_ACR_size_",size_summary_filt," n=",nrow(c_to_plot), sep=""))+ xlim(0,2000)


dev.off()

GRanges(seqnames =ACRs_classified$ACR_chr, ranges=IRanges(start=ACRs_classified$ACR_start,end=ACRs_classified$ACR_end),mcols=formcols)


cat_all <- GRanges(seqnames=Cat_all_acrs$V1, ranges=IRanges(start=Cat_all_acrs$V2,end=Cat_all_acrs$V3))

overlaps <- findOverlaps(cat_all, cat_all, type = "any")

# Calculate the proportion of each region that overlaps another region
overlap_proportions <- overlaps$width / lengths(cat_all)

# Plot the distribution of the proportions of overlap
ggplot(cat_all, aes(x = overlap_proportions)) +
  geom_histogram(bins = 20) +
  labs(x = "Proportion of region overlap")


#make column denoting whether an collapsed ACR overlaps atleast one ACR in a specific genotype, 0 means no , 1 means yes

Merged_ACRs.counts <- c
Merged_ACRs.counts$Overlap <- ifelse(Merged_ACRs.counts$V5 >= 1, 1, 0)



Overlap_summary <- Merged_ACRs.counts %>% group_by(V1,V2,V3) %>% summarize(Freq = sum(Overlap))

#acrs with freq of zero are due to initial filtering - ignore
Overlap_summary <- Overlap_summary %>% filter(Freq > 0)

#This produces table of unique ACRs and the genotypes they belong to
solos <- Overlap_summary %>% filter(Freq == 1)
solos_counts <- merge(Merged_ACRs.counts,solos , by=c("V1","V2","V3")) %>% filter(V5==1)

solos_table <-as.data.frame(table(solos_counts$V4))

z <- solos_counts %>% filter(V4 == "Interdonato_lemon")
Others <- solos_counts %>% filter(V4 != "Interdonato_lemon")

#cclem_genome <- read.delim("/bigdata/seymourlab/idiaz026/REFERENCES/Cclementina/Cclementine_genome.bed",header=F)[c(1:9),]

#custom.genome <- toGRanges(cclem_genome)
#kp <- plotKaryotype(genome = custom.genome)
#AS_ACRs <- toGRanges(z)
#all_acrs <- toGRanges(Others)

#kp <- plotKaryotype(genome=custom.genome, plot.type=2)

#kpPlotRegions(kp, data=AS_ACRs,col="red",data.panel=1)
#kpPlotRegions(kp, data=all_acrs,data.panel = 2)


pdf("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/2023-09-06_singleton_ACRs_bygeno.pdf")
ggplot(solos_table, aes(x=reorder(Var1, -Freq ), y=Freq )) + geom_bar(stat="identity") + coord_flip()+ theme_bw() + xlab("Genotype") + ylab("# Unique ACRs")

dev.off()


Unique_ACRs <- Overlap_summary %>% filter(Freq ==1)
Common_ACRs <- Overlap_summary %>% filter(Freq > 1 & Freq < 15)
High_freq_ACRs <- Overlap_summary %>% filter(Freq >= 16)

Unique_ACRs$Type <- "Unique"
Common_ACRs$Type <- "Common"
High_freq_ACRs$Type <-  "High_frequency"

ACRs_classified <- rbind(Unique_ACRs,Common_ACRs,High_freq_ACRs)[,c(-4)]

write.table(ACRs_classified, "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/ACRs_classified.bed" , sep ="\t", col.names=FALSE, row.names=FALSE,quote=FALSE)

write.table(Unique_ACRs[,c(1:3)], "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Unique_ACRs.bed", sep ="\t", col.names=FALSE, row.names=FALSE,quote=FALSE)

write.table(Common_ACRs[,c(1:3)], "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Common_ACRs.bed",sep="\t",col.names = FALSE, row.names=FALSE,quote =FALSE)

write.table(High_freq_ACRs[,c(1:3)], "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/High_frequency_ACRs.bed",sep="\t",col.names = FALSE, row.names=FALSE,quote =FALSE)

pdf("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACR_OVERLAP/2023-09-06_ACR_frequency.pdf")
ggplot(Overlap_summary,aes(x=Freq)) + geom_histogram(bins = 19,binwidth = 0.5) +  theme_classic() + scale_x_continuous(breaks = scales::pretty_breaks(n=19),limits=NULL)  + xlab("ACR frequency across 19 genotypes") + theme(text=element_text(size=15)) 

dev.off()

head(Overlap_summary)

#ACR present in one individual = unique, 2 - 16 = common, > 16 = high freq

#build matrix for upset plot
#have to use unfiltered data bc matrix dimensions
library(reshape2)

Merged_ACRs.counts <- read.delim("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACR_OVERLAP/Merged_ACRs.counts.bed", header=FALSE)

Merged_ACRs.counts$ACR_id <- paste(Merged_ACRs.counts$V1,Merged_ACRs.counts$V2,Merged_ACRs.counts$V3,sep="_")

Merged_ACRs.counts$V5 <- ifelse(Merged_ACRs.counts$V5 >= 1, 1, 0)

sub <- Merged_ACRs.counts[,c(4,5,6)]



ordered <- Merged_ACRs.counts[order(Merged_ACRs.counts$V4,Merged_ACRs.counts$ACR_id),]

names <- list(unique(ordered$ACR_id),unique(ordered$V4))

to_upset <- matrix(data=ordered$V5, ncol=19, nrow=length(unique(Merged_ACRs.counts$ACR_id)),byrow = F)

dimnames(to_upset) <- names

to_upset_df <- as.data.frame(to_upset)

upset(to_upset_df, sets = c("Microcitrus_australiasica", "Tango", "SCFS_citron", "Kao_panne_pummelo", "C_macrophylla"), sets.bar.color = "#56B4E9", group.by= "degree",order.by = "freq", empty.intersections = "on")


############################
##Focusing on select species 
############################

Merged_ACRs.counts <- read.delim("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACR_OVERLAP/Merged_ACRs.counts.species.bed", header=FALSE)

Cat_all_acrs <- read.delim("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACR_OVERLAP/Cat_species_acrs.sort.bed", header=FALSE)



Merged_ACRs.counts$Size <- Merged_ACRs.counts$V3 - Merged_ACRs.counts$V2


#Filter ACRs that are clustered - this leads to ACR boundaries expanding which biases ACR overlap

b <- Merged_ACRs.counts %>% filter(V5 > 0)

hist(b$V5, xlab = "# of ACRs in one sample overlapping consensus ACR")

b_to_plot <- unique(b[,c(1,2,3,6)],)

c <- Merged_ACRs.counts %>% filter(V5 < 2)

c_to_plot <- unique(c[,c(1,2,3,6)],)

pdf("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACR_OVERLAP/Species_level/2023-10-17_ACR.size.pdf")
#how does size of acrs change after collapsing overlapping ACRs.
Cat_all_acrs$Size <- Cat_all_acrs$V3 - Cat_all_acrs$V2
size_summary <- summary(Cat_all_acrs$Size)[3]
ggplot(Cat_all_acrs, aes(x=Size)) + geom_density() + theme_classic() + xlab("Size of all unmerged ACRs") + ggtitle(paste0("Median_ACR_size_",size_summary," n=",nrow(Cat_all_acrs), sep="")) + xlim(0,2000)


size_summary_merge <- summary(b_to_plot$Size)[3]
ggplot(b_to_plot, aes(x=Size)) + geom_density()  + theme_classic() + xlab("Size of all consensus ACRs") + ggtitle(paste0("Median_ACR_size_",size_summary_merge," n=",nrow(b_to_plot), sep="")) + xlim(0,2000)
summary(b_to_plot$Size)


num_filt <- dim(b_to_plot)[1] - dim(c_to_plot)[1]
hist(b$V5, xlab = "# of ACRs in one sample overlapping consensus ACR",main =paste0(num_filt," consensus ACRs filtered for being clustered",sep="" ))

size_summary_filt <- summary(c$Size)[3]
ggplot(c_to_plot, aes(x=Size)) + geom_density()  + theme_classic() + xlab("Size of all consensus, filtered ACRs") + ggtitle(paste0("Median_ACR_size_",size_summary_filt," n=",nrow(c_to_plot), sep=""))+ xlim(0,2000)


dev.off()


#make column denoting whether an collapsed ACR overlaps atleast one ACR in a specific genotype, 0 means no , 1 means yes

Merged_ACRs.counts <- c
Merged_ACRs.counts$Overlap <- ifelse(Merged_ACRs.counts$V5 >= 1, 1, 0)



Overlap_summary <- Merged_ACRs.counts %>% group_by(V1,V2,V3) %>% summarize(Freq = sum(Overlap))

#acrs with freq of zero are due to initial filtering - ignore
Overlap_summary <- Overlap_summary %>% filter(Freq > 0)

#This produces table of unique ACRs and the genotypes they belong to
solos <- Overlap_summary %>% filter(Freq == 1)
solos_counts <- merge(Merged_ACRs.counts,solos , by=c("V1","V2","V3")) %>% filter(V5==1)

solos_table <-as.data.frame(table(solos_counts$V4))

z <- solos_counts %>% filter(V4 == "Interdonato_lemon")
Others <- solos_counts %>% filter(V4 != "Interdonato_lemon")



uniquely_lost <-Overlap_summary %>% filter(Freq == 4)

uniquely_lost_counts <- merge(Merged_ACRs.counts,uniquely_lost, by=c("V1","V2","V3")) %>% filter(V5==0)

uniquely_lost_table <-as.data.frame(table(uniquely_lost_counts$V4))



pdf("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACR_OVERLAP/Species_level/2023-10-17_singleton_ACRs_bygeno.pdf")
ggplot(solos_table, aes(x=reorder(Var1, -Freq ), y=Freq )) + geom_bar(stat="identity") + coord_flip()+ theme_bw() + xlab("Genotype") + ylab("# Unique ACRs")

ggplot(uniquely_lost_table, aes(x=reorder(Var1, -Freq ), y=Freq )) + geom_bar(stat="identity") + coord_flip()+ theme_bw() + xlab("Genotype") + ylab("# Uniquely absent ACRs")

dev.off()


Unique_ACRs <- Overlap_summary %>% filter(Freq ==1)
Common_ACRs <- Overlap_summary %>% filter(Freq > 1 & Freq < 5)
High_freq_ACRs <- Overlap_summary %>% filter(Freq =5)

Unique_ACRs$Type <- "Unique"
Common_ACRs$Type <- "Common"
High_freq_ACRs$Type <-  "High_frequency"

ACRs_classified <- rbind(Unique_ACRs,Common_ACRs,High_freq_ACRs)[,c(-4)]

write.table(ACRs_classified, "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/ACRs_classified.bed" , sep ="\t", col.names=FALSE, row.names=FALSE,quote=FALSE)

write.table(Unique_ACRs[,c(1:3)], "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/Unique_ACRs.bed", sep ="\t", col.names=FALSE, row.names=FALSE,quote=FALSE)

write.table(Common_ACRs[,c(1:3)], "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/Common_ACRs.bed",sep="\t",col.names = FALSE, row.names=FALSE,quote =FALSE)

write.table(High_freq_ACRs[,c(1:3)], "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/High_frequency_ACRs.bed",sep="\t",col.names = FALSE, row.names=FALSE,quote =FALSE)

pdf("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACR_OVERLAP/Species_level/2023-10-17_ACR_frequency.pdf")
ggplot(Overlap_summary,aes(x=Freq)) + geom_histogram(bins = 5,binwidth = 0.5) +  theme_classic() + scale_x_continuous(breaks = scales::pretty_breaks(n=5),limits=NULL)  + xlab("ACR frequency across species") + theme(text=element_text(size=15)) 






library(reshape2)

Merged_ACRs.counts <- read.delim("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACR_OVERLAP/Merged_ACRs.counts.species.bed", header=FALSE)

Merged_ACRs.counts$ACR_id <- paste(Merged_ACRs.counts$V1,Merged_ACRs.counts$V2,Merged_ACRs.counts$V3,sep="_")

Merged_ACRs.counts$V5 <- ifelse(Merged_ACRs.counts$V5 >= 1, 1, 0)

sub <- Merged_ACRs.counts[,c(4,5,6)]



ordered <- Merged_ACRs.counts[order(Merged_ACRs.counts$V4,Merged_ACRs.counts$ACR_id),]

names <- list(unique(ordered$ACR_id),unique(ordered$V4))

to_upset <- matrix(data=ordered$V5, ncol=5, nrow=length(unique(Merged_ACRs.counts$ACR_id)),byrow = F)

dimnames(to_upset) <- names

to_upset_df <- as.data.frame(to_upset)

upset(to_upset_df, sets = c("Microcitrus_australiasica", "Seedless_kishu", "SCFS_citron", "Kao_panne_pummelo", "C_macrophylla"), sets.bar.color = "#56B4E9", group.by= "degree",order.by = "freq", empty.intersections = "on")


#plot uniquely lost ACRs 

dev.off()



