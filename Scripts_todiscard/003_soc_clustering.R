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

obj <- regModel(obj, verbose = T, nthreads=4, variates = "~genotype")
saveRDS(obj, file="/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Clustering/2023-09-06_GenoFilt_BEFORE_CLUST.rds")



print("reduce dims SVD")

obj <- reduceDims(obj, n.pcs=50, cor.max=0.7, verbose=T)


print("PROJECT UMAP")
obj <- projectUMAP(obj, verbose=T)



pdf("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Clustering/2023-09-14_UMAP.pdf")
par(mar=c(3,3,1,1))
plot(obj$UMAP, pch=16, cex=0.2, main="Socrates")
dev.off()

#print("CALL CLUSTERS")
obj <- callClusters(obj, res=0.6, verbose=T, k.near = 50, e.thresh =1.5)

pdf("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Clustering/2023-09-14_strict_umap_clusters.pdf")
plotUMAP(obj)
dev.off()


write.table(obj$Clusters, "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Clustering/2023-09-14_clusters.txt",col.names=F, row.names=F, quote=F,sep="\t")

write.table(obj$counts, "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Clustering/2023-09-14_counts.mtx", col.names=F, sep="\t",quote=F)


head(obj$UMAP)
b <- as.data.frame(obj$UMAP)



write.table(b,"/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Clustering/2023-09-14_GenoFilt_strict_umap_loadings.txt",sep="\t",quote=F,col.names=F)

b <- read.delim("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Clustering/2023-09-14_GenoFilt_strict_umap_loadings.txt")

obj <- coAccess(obj,genome=chr,byGroup=T,nthreads=2,verbose=T)

write.table(obj$coACRs, "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Clustering/2023-09-14_coACRs.txt",sep="\t", quote=F,row.names=T, col.names=T)

saveRDS(obj, file="/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Clustering/2023-09-14_GenoFilt_strict_Socrates_object.rds")

all <- read.table("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Clustering/2023-09-14_clusters.txt")


pdf("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Clustering/2023-09-06_umap_genotype_comp.pdf")

colnames(all) <- c("cellID","total","tss","acrs","ptmt","nSites","log10nSites","Genotype","UMAP1","UMAP2","clusterID")

genotypes <- unique(all$Genotype)


for (i in 1:length(genotypes)){
  write.table(all %>% filter(Genotype ==genotypes[i]) %>% select(cellID,Genotype), file=paste("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/SintoGenotypeBams/scATACFragments_files/",genotypes[i], ".postClust.barcodes.txt", sep = "" ), col.names=F , row.names = F, sep="\t", quote=F)
}

head(all)

ggplot(all, aes(x=UMAP1, y = UMAP2)) + geom_point(size=0.3,aes(color=log10nSites))+ theme_bw() + theme(legend.text=element_text(size=12), legend.key.height=unit(14, 'pt'))

ggplot(all, aes(x=UMAP1, y = UMAP2)) + geom_point(size=0.5,aes(color=Genotype))+ theme_bw() +scale_color_chameleon() + theme(legend.text=element_text(size=12), legend.key.height=unit(14, 'pt'))

dev.off()

clusters <- unique(all$clusterID)

clust_list <- list()
numcells <- c()

for (i in 1:length(clusters)) {
  clust_list[[i]] <- all %>% filter(clusterID == clusters[i])
  numcells[i] <- nrow(clust_list[[i]])

  }

avg_cells_per_cluster <- floor(mean(numcells))

pseudo1 <- c()
pseudo2 <- c()
pseudo3 <- c()

control1 <- c()
control2 <- c()
control3 <- c()

per_clust_sample <- floor(avg_cells_per_cluster / (length(clusters)-1))

#to srote %in%

pseudo1_store <- list()
pseudo2_store <- list()
pseudo3_store <- list()

control1_store <- list()
control2_store <- list()
control3_store <- list()

control_split <- list()


print("make pseudo reps and reference cell set")
for (i in 1:length(clusters)) {
  if (nrow(all %>% filter(clusterID == clusters[i])) >= avg_cells_per_cluster) {
    pseudo1[[i]] <- sample_n(all, avg_cells_per_cluster, replace = F)
    pseudo2[[i]] <- sample_n(all, avg_cells_per_cluster, replace = F)
    pseudo3[[i]] <- sample_n(all, avg_cells_per_cluster, replace=F)
  } else {
    pseudo1[[i]] <- sample_n(all, avg_cells_per_cluster, replace = T)
    pseudo2[[i]] <- sample_n(all, avg_cells_per_cluster, replace = T)
    pseudo3[[i]] <- sample_n(all, avg_cells_per_cluster, replace=T)
  }
  control1[[i]] <- sample_n((all %>% filter(clusterID == i)), per_clust_sample, replace=F)
  control2[[i]] <- sample_n((all %>% filter(clusterID == i)), per_clust_sample, replace=F)
  control3[[i]] <- sample_n((all %>% filter(clusterID == i)), per_clust_sample, replace=F)

  pseudo1[[i]]$rep <- paste("Cluster",clusters[i],"rep1",sep="_")
  pseudo2[[i]]$rep <- paste("Cluster",clusters[i],"rep2",sep="_")
  pseudo3[[i]]$rep <- paste("Cluster",clusters[i],"rep3",sep="_")

  pseudo1[[1]]$total_reads <- sum(pseudo1[[i]]$total)
  pseudo2[[1]]$total_reads <- sum(pseudo2[[i]]$total)
  pseudo3[[1]]$total_reads <- sum(pseudo3[[i]]$total)
}

control_1_df <- do.call(rbind,control1)
control_2_df <- do.call(rbind,control2)
control_3_df <- do.call(rbind,control3)

control_1_df$rep <- "Control_rep1"
control_2_df$rep <- "Control_rep2"
control_3_df$rep <- "Control_rep3"

control_1_df$total_reads <- sum(control_1_df$total)
control_2_df$total_reads <- sum(control_2_df$total)
control_3_df$total_reads <- sum(control_3_df$total)


cluster_spec_cont1 <- list()
cluster_spec_cont2 <- list()
cluster_spec_cont3 <- list()


print("Filtering controls so they dont include cells from cluster of interest")
for (i in 1:length(clusters)) {
  cluster_spec_cont1[[i]] <- control_1_df %>% filter(clusterID != i)
  cluster_spec_cont1[[i]]$rep <- paste("Control_for_Cluster_",i,"rep1",sep="")
  cluster_spec_cont2[[i]] <- control_2_df %>% filter(clusterID != i)
  cluster_spec_cont2[[i]]$rep <- paste("Control_for_Cluster_",i,"rep2",sep="")
  
  cluster_spec_cont3[[i]] <- control_3_df %>% filter(clusterID != i)
  
  cluster_spec_cont3[[i]]$rep <- paste("Control_for_Cluster_",i,"rep3",sep="")
  
}





for (i in 1:length(clusters)) {
  write.table(pseudo1[[i]], file=paste("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Counts_Gene_Activity/","Cluster",i,"_Rep1.metadata.txt",sep=""), sep = "\t",col.names=F, row.names=F,quote=F)

  write.table(pseudo2[[i]], file=paste("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Counts_Gene_Activity/","Cluster",i,"_Rep2.metadata.txt",sep=""), sep = "\t",col.names=F, row.names=F,quote=F)

  write.table(pseudo3[[i]], file=paste("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Counts_Gene_Activity/","Cluster",i,"_Rep3.metadata.txt",sep=""), sep = "\t",col.names=F, row.names=F,quote=F)

}

for (i in 1:length(clusters)) {
  write.table(cluster_spec_cont1[[i]], file=paste("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Counts_Gene_Activity/","Control_for_Cluster",i,"_Rep1.metadata.txt",sep=""), sep = "\t",col.names=F, row.names=F,quote=F)
  
  write.table(cluster_spec_cont2[[i]], file=paste("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Counts_Gene_Activity/","Control_for_Cluster",i,"_Rep2.metadata.txt",sep=""), sep = "\t",col.names=F, row.names=F,quote=F)
  
  write.table(cluster_spec_cont3[[i]], file=paste("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Counts_Gene_Activity/","Control_for_Cluster",i,"_Rep3.metadata.txt",sep=""), sep = "\t",col.names=F, row.names=F,quote=F)
  
}


#write.table(control_1_df, "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Counts_Gene_Activity/Control_rep1.metadata.txt", sep = "\t",col.names=F, row.names=F,quote=F)
#write.table(control_2_df, "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Counts_Gene_Activity/Control_rep2.metadata.txt", sep = "\t",col.names=F, row.names=F,quote=F)
#write.table(control_3_df, "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Counts_Gene_Activity/Control_rep3.metadata.txt", sep = "\t",col.names=F, row.names=F,quote=F)


head(control_1_df)


# #TESTING


# obj1 <- callClusters(obj, res=0.6, verbose=T)

# obj2 <- callClusters(obj, res=0.6, verbose=T, k.near = 30, e.thresh =1)

# obj3 <- callClusters(obj, res=0.6, verbose=T, k.near = 30, e.thresh =1.5)

# obj4 <- callClusters(obj, res=0.6, verbose=T, k.near = 30, e.thresh =0.5)

# pdf("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Clustering/clusterz-test_res6.pdf")

# plotUMAP(obj)
# plotUMAP(obj1)
# plotUMAP(obj2)
# plotUMAP(obj3)
# plotUMAP(obj4)
# dev.off()

# obj1 <- callClusters(obj, res=0.8, verbose=T)

# obj2 <- callClusters(obj, res=1.0, verbose=T)

# obj3 <- callClusters(obj, res=1.2, verbose=T)

# obj4 <- callClusters(obj, res=1.6, verbose=T)

# obj5 <- callClusters(obj, res=1.8, verbose=T)
# obj6 <- callClusters(obj, res=2.0, verbose=T)

# pdf("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Clustering/clusterz-test_res_test.pdf")


# plotUMAP(obj1)
# plotUMAP(obj2)
# plotUMAP(obj3)
# plotUMAP(obj4)
# plotUMAP(obj5)
# plotUMAP(obj6)

# dev.off()

# obj1 <- callClusters(obj, res=0.6, verbose=T, k.near = 20, e.thresh =1.5)

# obj2 <- callClusters(obj, res=0.6, verbose=T, k.near = 30, e.thresh =1.5)
# obj3 <- callClusters(obj, res=0.6, verbose=T, k.near = 40, e.thresh =1.5)
# obj4 <- callClusters(obj, res=0.6, verbose=T, k.near = 50, e.thresh =1.5)
# obj5 <- callClusters(obj, res=0.6, verbose=T, k.near = 60, e.thresh =1.5)


# pdf("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Clustering/clusterz-test_knear_test.pdf")

# plotUMAP(obj1)
# plotUMAP(obj2)
# plotUMAP(obj3)
# plotUMAP(obj4)
# plotUMAP(obj5)

# dev.off()

# obj1 <- callClusters(obj, res=0.6, verbose=T, k.near = 50, e.thresh =1.5,cl.method=1)
# obj2 <- callClusters(obj, res=0.6, verbose=T, k.near = 50, e.thresh =1.5,cl.method=2)
# obj3 <- callClusters(obj, res=0.6, verbose=T, k.near = 50, e.thresh =1.5,cl.method=3)
# obj4 <- callClusters(obj, res=0.6, verbose=T, k.near = 50, e.thresh =1.5,cl.method=4)


# pdf("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Clustering/clusterz-test_method_test.pdf")

# plotUMAP(obj1)
# plotUMAP(obj2)
# plotUMAP(obj3)

# dev.off()
