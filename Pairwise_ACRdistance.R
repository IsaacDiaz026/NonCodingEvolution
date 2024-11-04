### 2022-05-08 
#Calculate pairwise differences in ACRs
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(tidyr)
library(ggforce)
library(dplyr)

setwd("/bigdata/seymourlab/idiaz026/NonCodingEvolution/Plink/")

ALL_ACRs.sort <- read.delim("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/ACRs_classified.sort.bed", header=FALSE)

ALL_ACRs.sort$size <- ALL_ACRs.sort$V3 - ALL_ACRs.sort$V2

TOTAL_ACR <- sum(ALL_ACRs.sort$size)


blues <- colorRampPalette(c('dark blue', 'light blue'))
greens <- colorRampPalette(c('dark green', 'light green'))
reds <- colorRampPalette(c('pink', 'dark red'))

Pairwise_dist <- read.delim("VCF_75_miss.diff.sdiff.summary")


scATAC_stats <- read.table("/bigdata/seymourlab/idiaz026/Data/Floral_scATAC/scATAC_stats_open_chrom.txt", header = T, sep = "\t")

scATAC_stats <- scATAC_stats %>% filter(genotype != "Interdonato_lemon")

scATAC_stats$readspercell <- scATAC_stats$Number.of.Reads/scATAC_stats$Number.of.Nuclei
600000 / TOTAL_ACR

diff_table <- list()
diff_table2 <- list()
diff_matrix <- list()

vcf_type <- c("75_miss")

Difference<- noquote(as.numeric(rep(0, 19)))
X.IID1 <- scATAC_stats$Sample_ID
IID2 <- scATAC_stats$Sample_ID
add_fake_cor <- as.data.frame(cbind(X.IID1,IID2,Difference))


diff_table <- Pairwise_dist
diff_table$Difference <- (diff_table$DIFF_CT)
diff_table <- diff_table[,c(1,2,5)]
#add fake 1 to 1 corr because plink doesnt output this
diff_table2 <- diff_table[,c(2,1,3)]
colnames(diff_table2) <- c("X.IID1","IID2","Difference")
 
diff_table <- as.data.frame(rbind(diff_table,diff_table2))
diff_table <- as.data.frame(rbind(diff_table, add_fake_cor))
diff_table$Difference <- as.numeric(diff_table$Difference)
 
diff_table <- spread(diff_table,IID2,Difference)[,-1]
 
diff_matrix <- as.matrix(diff_table)
row.names(diff_matrix) <- scATAC_stats[,c(1)]

colnames(diff_matrix) <- scATAC_stats$genotype
 


hclust_rows <- as.dendrogram(hclust(dist(diff_matrix)))
hclust_cols <- as.dendrogram(hclust(dist(t(diff_matrix))))


plot(hclust_cols)

pdf("/bigdata/seymourlab/idiaz026/Plots/Floral_scATAC/Summarize_ATAC/ACR_pairwise_dist/Pairwise_dist.heatmap.pdf")
heatmap.2(diff_matrix, col = brewer.pal(9,"Greens"), margins = c(10, 10), density.info = "none", lhei=c(2, 8), trace= "none",dendrogram = "row", Rowv = hclust_rows, Colv=hclust_cols, main=vcf_type[1] )

dev.off()

n <- length(diff_matrix[,1])
diff_df <- vector()


for (i in 1:19) {
 diff_df[i] <- sort(diff_matrix[,i],partial = n-17)[n-17]
 
}

write.table(closest_w_meta, "/bigdata/seymourlab/idiaz026/Plots/Floral_scATAC/Summarize_ATAC/ACR_pairwise_dist/most_similar_genetic_distance.txt",col.names=T, sep="\t",row.names=T)

closest_w_meta <- cbind(scATAC_stats,diff_df)


ggplot(closest_w_meta, aes(x=reorder(genotype, -diff_df),y=diff_df,fill=Market_Type)) + geom_bar(stat="identity") + coord_flip()
barplot(closest_w_meta$diff_df)




#for binom distribution we need probability of success which is equal to $Diff column, 
#we need number of attempts = reads per cell * FRIP * read length 
#and we need the number of desired successes we will do a range of 

FRIPS <- read.table("bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-05-05_Downsample_ACR_analysis/FRIP/Collate_peaks_allFrip.txt", header=T)

FRIPS$file <- gsub("\\/(.*?).sort.bam","\\1",FRIPS$file)

FRIPS <- FRIPS[,c(1,3)]

Full_metadata <- merge(closest_w_meta, FRIPS, by.x =c("genotype"), by.y=c("file"))

Full_metadata$Number_attempts <- Full_metadata$readspercell * (Full_metadata$percent/100)

View(diff_matrix[[1]])


#make columns for probabillity of sampling varying numbers of SNPs
pbinom(q=1000,size= floor(Full_metadata$Number_attempts), prob= Full_metadata$diff_df, lower.tail = F)

probs <- list()
sub <- list()
x <- c(1:1500)

for (i in (1:20)) {
 sub[[i]] <- Full_metadata %>% filter(genotype == Full_metadata[i,1])
 probs[[i]] <- data.frame(pbinom(x,size=floor(sub[[i]]$Number_attempts),prob=sub[[i]]$diff_df), lower.tail=T)
 probs[[i]]$genotype <- sub[[i]]$genotype
 probs[[i]]$k <- c(1:1500)
 colnames(probs[[i]]) <- c("Probabillity", "True","Genotype", "Number_of_Successes")
}

ggplot(probs[[9]], aes(x= Number_of_Successes, y = Probabillity)) + geom_point() + theme_classic()

head(probs[[11]])

b <- (1:5000)
prob <- data.frame(pbinom(b,size=floor(sub[[9]]$Number_attempts),prob=sub[[9]]$diff_df), lower.tail=T)
prob$k <- c(1:5000)
colnames(prob) <- c("prob","true","k")

ggplot(prob, aes(x=k, y =prob ))+ geom_point()

All_binom <- do.call(rbind, probs)

ggplot(All_binom, aes(x=Number_of_Successes, y = Probabillity, color = Genotype)) + geom_point() + theme_classic() + xlim(0,900)

All_binom_sub <- All_binom %>% filter(Genotype == "Tango" | Genotype == "Koster_mandarin" | Genotype == "Interdonato_lemon" | Genotype == "Washington_navel" | Genotype == "Washington_navel")

png("/bigdata/seymourlab/idiaz026/Plots/Floral_scATAC/Summarize_ATAC/ACR_pairwise_dist/ACR_binom.sub.png")
ggplot(All_binom_sub, aes(x=Number_of_Successes, y = Probabillity, color = Genotype)) + geom_point() + theme_classic() + xlim(0,50) + theme(text=element_text(size=20))

dev.off()

Full_metadata$Confidence <- Full_metadata$Number_attempts * Full_metadata$diff_df

png("/bigdata/seymourlab/idiaz026/Plots/Floral_scATAC/Summarize_ATAC/ACR_pairwise_dist/Nuclei_classif_confidence.png")

ggplot(Full_metadata, aes(x=reorder(genotype,-Confidence), y = (Confidence), fill =Market_Type)) + geom_bar(stat="identity") + coord_flip() + theme_classic() + ylab("Confidence") + xlab("Genotype") + theme(text=element_text(size=20)) + scale_fill_brewer(palette="Set3")

dev.off()

b <- melt(Full_metadata[,c(1,8,9,10)])
ggplot(b, aes(fill=variable, y = value, x=genotype))+ geom_bar(position="dodge",stat="identity") +coord_flip() + theme_classic() + theme(text=element_text(size=20)) + scale_fill_brewer(palette="Set2") + ylab("Informative SNPs / Open chromatin") + xlab("Genotype")



