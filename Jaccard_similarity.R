#Plot heatmap of jaccard similarity between genotype
#2022-02-21
setwd("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_by_geno/")

library(ggplot2)
library(RColorBrewer)
library(gplots)
library(tidyr)
library(ggforce)

#setup colors
blues <- colorRampPalette(c('dark blue', 'light blue'))
greens <- colorRampPalette(c('dark green', 'light green'))
reds <- colorRampPalette(c('pink', 'dark red'))

pairwise.ATAC.shortnames <- read.delim("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_by_geno/Downsampled_genrich/pairwise.ATAC.shortnames.txt", header=FALSE)


jaccard_table <-spread(pairwise.ATAC.shortnames,V2,V3)[,-1]

jaccard_matrix <- as.matrix(jaccard_table)
row.names(jaccard_matrix) <- colnames(jaccard_table)

species <- c("Seedless_kishu_sort.CLEAN.bed","C_macrophylla.sort.CLEAN.bed",
             "Kao_panne_pummelo.sort.CLEAN.bed","Microcitrus_australiasica.sort.CLEAN.bed",
             "SCFS_citron.sort.CLEAN.bed")

species_jaccard <- jaccard_matrix[(rownames(jaccard_matrix) %in% species),]
species_jaccard <- species_jaccard[,(colnames(species_jaccard) %in% species)]

pdf("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_by_geno/2023-09-05_Jaccard_similarity_ACRs_orig.pdf")
heatmap.2(jaccard_matrix, col = brewer.pal(9,"Blues"), margins = c(10, 10), density.info = "none", lhei=c(2, 8), trace= "none")
heatmap.2(species_jaccard, col = brewer.pal(9,"Blues"), margins = c(10, 10), density.info = "none", lhei=c(2, 8), trace= "none")

dev.off()

