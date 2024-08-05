#!/usr/bin/env Rscript

library(devtools)

library(Socrates)
library(dplyr)
library(doSNOW)
library(ggplot2)
library(chameleon)

geno <- read.table("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/SintoGenotypeBams/2023-09-06_all_snps_barcodesClean_minusInter.txt")

colnames(geno) <-c("cellID","genotype")

meta <- read.table("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/2023-08-24_cell_meta_beforeisCell.txt",header=T)

print("merging meta with genotype data")
meta <- merge(meta, geno , by=c("cellID"))

print(head(meta))

row.names(meta) <- meta$cellID
write.table(meta,"/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Clustering/2023-09-06_Old_meta_minusInter.txt", col.names =T, row.names=T, quote=F, sep="\t")


genes <- read.table("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/WN_HAPA_genes_modCTG.bed")

colnames(genes) <- c("chr","start","end","strand","geneid")

pos <- subset(genes, strand == "+")
pos$end <- pos$end + 100
pos$start <- pos$start - 500


neg <- subset(genes, strand == "-")
neg$start <- neg$start - 100
neg$end <- neg$end + 500

genes <- rbind(pos,neg)

write.table(genes, "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Counts_Gene_Activity/Genes_upAndDownstream.bed", sep ="\t", col.names=F, row.names=F, quote=F)


meta <- "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Clustering/2023-09-06_Old_meta_minusInter.txt"

print("load matrix from combined matrices filtered by true cells")

matrices <- "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/All_ATAC_peaks.sparse.gz"

chr <- "/bigdata/seymourlab/idiaz026/REFERENCES/PWN_HAPA/WN_HAPA_modContig_genomeFile.txt"

print("create socrates obj")
obj <- loadSparseData(input=matrices, meta=meta, verbose=T)


print("estimate log10 number of accessible regions per cell")
cell.counts <- Matrix::colSums(obj$counts)

print("estimate peak accessibility frequency across cells")
site.freq <- Matrix::rowMeans(obj$counts)

print("plot distributions")

pdf("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Clustering/2023-09-06_acrsPerCell.pdf")

layout(matrix(c(1:2), ncol=2))
par(mar=c(3,3,1,1))
plot(density(cell.counts), main="log10 cell counts", log="x")


abline(v=1000, col="red")
plot(density(site.freq), main="peak accessibility frequency", log="x")

dev.off()

obj <- cleanData(obj, verbose=T, min.t=0.03, min.c=1000)

print(head(obj$meta))
print("regress out geno")

obj <- regModel(obj, verbose = T, nthreads=4)
saveRDS(obj, file="/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Clustering/2023-09-17_NOTGENOTYPE_CORRECTED_BEFORE_CLUST.rds")
