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
library(reshape2)




pairwise.ATAC.shortnames <- fread(args[1])
genrich <- fread(args[2])

Merged_ACRs.counts <- fread(args[3]) # Consensus ACRs
Cat_all_acrs <- fread(args[4]) #All ACRs

PREFIX <- args[5] # specify "Total" for all samples, "Species_level" for species


#setup colors
blues <- colorRampPalette(c('dark blue', 'light blue'))
greens <- colorRampPalette(c('dark green', 'light green'))
reds <- colorRampPalette(c('pink', 'dark red'))

jaccard_table <-spread(pairwise.ATAC.shortnames,V2,V3)[,-1]

jaccard_matrix <- as.matrix(jaccard_table)
row.names(jaccard_matrix) <- colnames(jaccard_table)

#make matrix with only species
species <- c("Algerian_clementine.sort.CLEAN.bed","C_macrophylla.sort.CLEAN.bed",
             "Kao_panne_pummelo.sort.CLEAN.bed","Microcitrus_australiasica.sort.CLEAN.bed",
             "SCFS_citron.sort.CLEAN.bed")

species_jaccard <- jaccard_matrix[(rownames(jaccard_matrix) %in% species),]
species_jaccard <- species_jaccard[,(colnames(species_jaccard) %in% species)]

pdf("Jaccard_similarity_ACRs.pdf")
heatmap.2(jaccard_matrix, col = brewer.pal(9,"Blues"), margins = c(10, 10), density.info = "none", lhei=c(2, 8), trace= "none")
heatmap.2(species_jaccard, col = brewer.pal(9,"Blues"), margins = c(10, 10), density.info = "none", lhei=c(2, 8), trace= "none")

dev.off()


pdf("ACR.pos.overlap.pdf")

plot(genrich[,1], genrich_loose[,2] / sum(genrich_loose[,2]), 'h',
     col="darkred",lwd=4, xlab="No. of assayed genotypes",
     ylab="Fraction of bases",main="Genrich_0.05")

dev.off()

Merged_ACRs.counts <- fread(args[3])

Cat_all_acrs <- fread(args[4])

Merged_ACRs.counts$Size <- Merged_ACRs.counts$V3 - Merged_ACRs.counts$V2


#Filter ACRs that are clustered - this leads to ACR boundaries expanding which biases ACR overlap
b <- Merged_ACRs.counts %>% filter(V5 > 0)
b_to_plot <- unique(b[,c(1,2,3,6)],)

#filter clustered 
c <- Merged_ACRs.counts %>% filter(V5 < 2)
c_to_plot <- unique(c[,c(1,2,3,6)],)


pdf("ACR.size.pdf")
#how does size of acrs change after collapsing overlapping ACRs.
hist(b$V5, xlab = "# of ACRs in one sample overlapping consensus ACR")

Cat_all_acrs$Size <- Cat_all_acrs$V3 - Cat_all_acrs$V2
size_summary <- summary(Cat_all_acrs$Size)[3]
ggplot(Cat_all_acrs, aes(x=Size)) + geom_density() + theme_classic()
+ xlab("Size of all unmerged ACRs") + 
ggtitle(paste0("Median_ACR_size_",size_summary," n=",nrow(Cat_all_acrs), sep="")) + xlim(0,2000)


size_summary_merge <- summary(b_to_plot$Size)[3]
ggplot(b_to_plot, aes(x=Size)) + geom_density()  + theme_classic() + 
xlab("Size of all consensus ACRsbefore removing clustered") + 
ggtitle(paste0("Median_ACR_size_",size_summary_merge," n=",nrow(b_to_plot), sep=""))+ xlim(0,2000)

summary(b_to_plot$Size)


num_filt <- dim(b_to_plot)[1] - dim(c_to_plot)[1]
hist(b$V5, xlab = "# of ACRs in one sample overlapping consensus ACR",main =paste0(num_filt," consensus ACRs filtered for being clustered",sep="" ))

size_summary_filt <- summary(c_to_plot$Size)[3]
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


#####

#percent overlap 
cat_all <- GRanges(seqnames=Cat_all_acrs$V1, ranges=IRanges(start=Cat_all_acrs$V2,end=Cat_all_acrs$V3))

cat_all <- setNames(cat_all,Cat_all_acrs$ID)

consensus <- GRanges(seqnames=Overlap_summary$V1, ranges=IRanges(start=Overlap_summary$V2,end=Overlap_summary$V3))

values(consensus) <- DataFrame(Freq=Overlap_summary$Freq)

consensus <- setNames(consensus,paste(Overlap_summary$V1,Overlap_summary$V2,Overlap_summary$V3,Overlap_summary$Freq,sep="_"))


hits <- findOverlaps(cat_all, consensus, minoverlap = 1)

overlaps <- pintersect(cat_all[queryHits(hits)], consensus[subjectHits(hits)])

percentOverlap <- (width(overlaps) / width(consensus[subjectHits(hits)])) * 100

median_overlap <- summary(percentOverlap)[3]
# Calculate the proportion of each region that overlaps another region

# Plot the distribution of the proportions of overlap
hist(percentOverlap, main=paste("Overlap all ACRs to consensus ACR","median",median_overlap))


#get freq info for all acrs from this overlap

ranges <- subsetByOverlaps(cat_all,consensus)
hits <- findOverlaps(cat_all, consensus)

rsid <- CharacterList(split(names(consensus)[subjectHits(hits)],
                              queryHits(hits)))

freq <- CharacterList(split(consensus$Freq[subjectHits(hits)],
                                  queryHits(hits)))

mcols(ranges) <- DataFrame(mcols(ranges), rsid, freq)


names(ranges) <- NULL
df <- unique(as.data.frame(ranges)[,c(1,2,3,4,7)])

unique_all <- df %>% filter(freq == 1)
common_all <- df %>% filter(freq > 1 & freq < 15)
high_all <- df %>% filter(freq > 15)

unique_all$Type <- "Unique"
common_all$Type <- "common"
high_all$Type <- "High"

cat_allfreq <- rbind(unique_all,common_all,high_all)

ggplot(cat_allfreq, aes(x=Type, y=log(width))) + geom_boxplot()

cat_allfreq$Size <- log(cat_allfreq$width)


pdf("singleton_ACRs_bygeno.pdf")
ggplot(solos_table, aes(x=reorder(Var1, -Freq ), y=Freq )) + geom_bar(stat="identity") + coord_flip()+ theme_bw() + xlab("Genotype") + ylab("# Unique ACRs")

dev.off()


Unique_ACRs <- Overlap_summary %>% filter(Freq ==1)
Common_ACRs <- Overlap_summary %>% filter(Freq > 1 & Freq < 15)
High_freq_ACRs <- Overlap_summary %>% filter(Freq >= 16)

Unique_ACRs$Type <- "Unique"
Common_ACRs$Type <- "Common"
High_freq_ACRs$Type <-  "High_frequency"

ACRs_classified <- rbind(Unique_ACRs,Common_ACRs,High_freq_ACRs)[,c(-4)]

# Export consensus ACRs 

write.table(ACRs_classified, paste("_ACRs_classified/",PREFIX,"_ACRs_classified.bed",sep="") , sep ="\t", col.names=FALSE, row.names=FALSE,quote=FALSE)

write.table(Unique_ACRs[,c(1:3)], paste("_ACRs_classified/",PREFIX,"Unique_ACRs.bed",sep=""), sep ="\t", col.names=FALSE, row.names=FALSE,quote=FALSE)

write.table(Common_ACRs[,c(1:3)], paste("_ACRs_classified/",PREFIX,"Common_ACRs.bed",sep=""),sep="\t",col.names = FALSE, row.names=FALSE,quote =FALSE)

write.table(High_freq_ACRs[,c(1:3)], paste("_ACRs_classified/",PREFIX,"High_frequency_ACRs.bed",sep=""),sep="\t",col.names = FALSE, row.names=FALSE,quote =FALSE)

#Export ACRs prior to consensus building "Cat_all" 
write.table(cat_allfreq[,c(1,2,3,6)], paste("_ACRs_classified/",PREFIX, "CAT_ALLACRs_classified.bed",sep="") , sep ="\t", col.names=FALSE, row.names=FALSE,quote=FALSE)

write.table(unique_all[,c(1:3)], paste("_ACRs_classified/",PREFIX,"CatALLUnique_ACRs.bed",sep=""), sep ="\t", col.names=FALSE, row.names=FALSE,quote=FALSE)

write.table(common_all[,c(1:3)], paste("_ACRs_classified/",PREFIX, "CatALLCommon_ACRs.bed",sep=""),sep="\t",col.names = FALSE, row.names=FALSE,quote =FALSE)

write.table(high_all[,c(1:3)], paste("_ACRs_classified/",PREFIX,"CatALL_High_frequency_ACRs.bed",sep=""),sep="\t",col.names = FALSE, row.names=FALSE,quote =FALSE)


pdf(paste("ACRs_classified/",PREFIX,"ACR_frequency.pdf",sep=""))
ggplot(Overlap_summary,aes(x=Freq)) + geom_histogram(bins = 19,binwidth = 0.5) + theme_classic() + scale_x_continuous(breaks = scales::pretty_breaks(n=19),limits=NULL)  + xlab("ACR frequency across 19 genotypes") + theme(text=element_text(size=15)) 

dev.off()

head(Overlap_summary)

#ACR present in one individual = unique, 2 - 16 = common, > 16 = high freq

#build matrix for upset plot
#have to use unfiltered data bc matrix dimensions

Merged_ACRs.counts <- fread(args[3])

n <- length(unique(Merged_ACRs.counts$V4))

Merged_ACRs.counts$ACR_id <- paste(Merged_ACRs.counts$V1,Merged_ACRs.counts$V2,Merged_ACRs.counts$V3,sep="_")

Merged_ACRs.counts$V5 <- ifelse(Merged_ACRs.counts$V5 >= 1, 1, 0)

sub <- Merged_ACRs.counts[,c(4,5,6)]


ordered <- Merged_ACRs.counts[order(Merged_ACRs.counts$V4,Merged_ACRs.counts$ACR_id),]

names <- list(unique(ordered$ACR_id),unique(ordered$V4))

to_upset <- matrix(data=ordered$V5, ncol=n, nrow=length(unique(Merged_ACRs.counts$ACR_id)),byrow = F)

dimnames(to_upset) <- names

to_upset_df <- as.data.frame(to_upset)

upset(to_upset_df, sets = c("Microcitrus_australiasica", "Algerian_clementine", "SCFS_citron", "Kao_panne_pummelo", "C_macrophylla"), sets.bar.color = "#56B4E9", group.by= "degree",order.by = "freq", empty.intersections = "on")






## for next script
uniquely_lost <-Overlap_summary %>% filter(Freq == 4)

uniquely_lost_counts <- merge(Merged_ACRs.counts,uniquely_lost, by=c("V1","V2","V3")) %>% filter(V5==0)

uniquely_lost_table <-as.data.frame(table(uniquely_lost_counts$V4))


