library(Socrates)
library(dplyr)



input <- "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/All_ATAC_peaks.sparse.gz"

metafile <- "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Clustering/2023-09-14_clusters.txt"

threads=5
a <- read.table(input)

zm <- read.table("/bigdata/seymourlab/idiaz026/REFERENCES/PWN_HAPA/WN_HAPA_modContig_genomeFile.txt")

genes <- read.table("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/WN_HAPA_genes_modCTG.bed")

colnames(genes) <-c("chr","start","end","strand","transcript")

meta <- read.table(metafile)

coaccess <- "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Clustering/2023-09-06_coACRs.txt"

colnames(meta) <- c("cellID","total","tss","acrs","ptmt","nSites","log10nSites","Genotype","umap1","umap2","clusterID")

row.names(meta) <- meta$cellID


out <- "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Cicero_Gene_Activity/"

library(cicero)
library(Matrix)
library(parallel)
library(doSNOW)
library(methods)
library(tcltk)
library(iterators)
library(itertools)


source("moreCiceroUTILS.R")
###############################################################################
## Cicero trajectories
###############################################################################


###################################################################################################
## Cluster cells 
###################################################################################################

# create cicero object
message(" ... Creating CDS")
a <- a[as.character(a$V2) %in% rownames(meta),]
#a$V1 <- droplevels(a$V1)
#a$V2 <- droplevels(a$V2)
shuf <- a
shuf$V1 <- shuf$V1[sample(length(shuf$V1))]
shuf$V2 <- shuf$V2[sample(length(shuf$V2))]
cds <- make_atac_cds(a, binarize=T)
shufcds <- make_atac_cds(shuf, binarize=T)

c.colSums <- Matrix::colSums(exprs(cds))
c.rowSums <- Matrix::rowSums(exprs(cds))
#s.colSums <- Matrix::colSums(exprs(shufcds))
#s.rowSums <- Matrix::rowSums(exprs(shufcds))

# check shuffled
#message("   * orig. matrix = ",nrow(cds), " | ", ncol(cds))
# message("   * shuf. matrix = ",nrow(shufcds), " | ", ncol(shufcds))
# message("   * orig. matrix colSums = ", paste(c.colSums[1:5], collapse=", "))
# message("   * shuf. matrix colSums = ", paste(s.colSums[1:5], collapse=", "))
# message("   * orig. matrix rowSums = ", paste(c.rowSums[1:5], collapse=", "))
# message("   * shuf. matrix rowSums = ", paste(s.rowSums[1:5], collapse=", "))

# add metadata, filter, and run TFIDF/library regression/batch effect removal
#cds <- cds[,colnames(cds) %in% rownames(meta)]
#shufcds <- shufcds[,colnames(cds) %in% rownames(meta)]
pData(cds) <- meta[colnames(exprs(cds)),]
#pData(shufcds) <- meta[colnames(exprs(shufcds)),]
cds <- cds[Matrix::rowSums(exprs(cds))>0,]
cds <- cds[,Matrix::colSums(exprs(cds))>0]
#shufcds <- shufcds[Matrix::rowSums(exprs(shufcds))>0,]
#shufcds <- shufcds[,Matrix::colSums(exprs(shufcds))>0]

# process basic
cds <- detectGenes(cds)
#shufcds <- detectGenes(shufcds)
cds <- estimateSizeFactors(cds)
#shufcds <- estimateSizeFactors(shufcds)



# .loadMeta <- function(cds, meta, clusterID){

#         # svd and raw
#         ids <- colnames(cds)
#         meta <- meta[ids,]
#         meta$UMAP1 <- meta$umap1
#         meta$UMAP2 <- meta$umap2
#         ids.2 <- rownames(meta)
#         cds <- cds[,colnames(cds) %in% ids.2]
#         umap.d <- t(meta[,c("UMAP1","UMAP2")])

#         # UMAP output
#         cds@reducedDimA <- umap.d
#         cds@reducedDimS <- umap.d
#         cds@dim_reduce_type <- "UMAP"
#         pData(cds)$Cluster <- as.factor(meta$clusterID)
#         # return data
#         return(cds)
#     }
 
# cds <- .loadMeta(cds,meta,1)


meta2 <- pData(cds)
meta2$Cluster <- as.character(meta2$clusterID)
print(table(meta2$Cluster))
clusts <- unique(meta2$Cluster)
cell_ids <- c()

# iterate
its <- 0

# foreach parameters
cl <- makeSOCKcluster(threads)
registerDoSNOW(cl)
tasks <- length(clusts)
pb <- txtProgressBar(max = tasks, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
package.labs <- c("cicero", "Matrix")




build_composite_gene_activity_matrix <- function(input_cds,
                                                 site_weights,
                                                 cicero_cons_info,
                                                 dist_thresh=250000,
                                                 coaccess_cutoff=0.25) {
    accessibility_mat <- exprs(input_cds)
    promoter_peak_table <- fData(input_cds)
    promoter_peak_table$peak <- as.character(row.names(promoter_peak_table))
    promoter_peak_table <-
        promoter_peak_table[!is.na(promoter_peak_table$gene),]
    promoter_peak_table <- promoter_peak_table[,c("peak", "gene")]
    promoter_peak_table$gene <- as.character(promoter_peak_table$gene)

    # Make site_weight matrix
    site_names <- names(site_weights)
    site_weights <- as(Matrix::Diagonal(x=as.numeric(site_weights)),
                      "sparseMatrix")
    row.names(site_weights) <- site_names
    colnames(site_weights) <- site_names

    # Find distance between cicero peaks. If distance already calculated, skip
    if ("dist" %in% colnames(cicero_cons_info) == FALSE) {
        Peak1_cols <- split_peak_names(cicero_cons_info$Peak1)
        Peak2_cols <- split_peak_names(cicero_cons_info$Peak2)
        Peak1_bp <- round((as.integer(Peak1_cols[,3]) +
                          as.integer(Peak1_cols[,2])) / 2)
        Peak2_bp <- round((as.integer(Peak2_cols[,3]) +
                          as.integer(Peak2_cols[,2])) / 2)
        cicero_cons_info$dist <- abs(Peak2_bp - Peak1_bp)
    }

    # Get connections between promoters and distal sites above coaccess
    # threshold
    nonneg_cons <-
        cicero_cons_info[(cicero_cons_info$Peak1 %in%
                          promoter_peak_table$peak |
                          cicero_cons_info$Peak2 %in%
                          promoter_peak_table$peak) &
                          cicero_cons_info$coaccess >= coaccess_cutoff &
                          cicero_cons_info$dist < dist_thresh,]
    nonneg_cons <- nonneg_cons[,c("Peak1", "Peak2", "coaccess")]
    nonneg_cons <- nonneg_cons[!duplicated(nonneg_cons),]

    nonneg_cons$Peak1 <- as.character(nonneg_cons$Peak1)
    nonneg_cons$Peak2 <- as.character(nonneg_cons$Peak2)

    nonneg_cons <- rbind(nonneg_cons,
                        data.frame(Peak1=unique(promoter_peak_table$peak),
                                   Peak2=unique(promoter_peak_table$peak),
                                   coaccess=0))

    # Make square matrix of connections from distal to proximal
    distal_connectivity_matrix <- make_sparse_matrix(nonneg_cons,
                                                    x.name="coaccess")

    # Make connectivity matrix of promoters versus all
    promoter_conn_matrix <-
        distal_connectivity_matrix[unique(promoter_peak_table$peak),]

    # Get list of promoter and distal sites in accessibility mat
    promoter_safe_sites <- intersect(rownames(promoter_conn_matrix),
                                     row.names(accessibility_mat))
    distal_safe_sites <- intersect(colnames(promoter_conn_matrix),
                                     row.names(accessibility_mat))
    distal_safe_sites <- setdiff(distal_safe_sites, promoter_safe_sites)

    # Get accessibility info for promoters
    promoter_access_mat_in_cicero_map <- accessibility_mat[promoter_safe_sites,, drop=FALSE]

    # Get accessibility for distal sites
    distal_activity_scores <- accessibility_mat[distal_safe_sites,, drop=FALSE]

    # Scale connectivity matrix by site_weights
    scaled_site_weights <- site_weights[distal_safe_sites,distal_safe_sites, drop=FALSE]
    total_linked_site_weights <- promoter_conn_matrix[,distal_safe_sites, drop=FALSE] %*%
        scaled_site_weights
    total_linked_site_weights <- 1/Matrix::rowSums(total_linked_site_weights,
                                                na.rm=TRUE)
    total_linked_site_weights[is.finite(total_linked_site_weights) == FALSE] <- 0
    total_linked_site_weights[is.na(total_linked_site_weights)] <- 0
    total_linked_site_weights[is.nan(total_linked_site_weights)] <- 0
    total_linked_site_weights <- Matrix::Diagonal(x=total_linked_site_weights)
    scaled_site_weights <- total_linked_site_weights %*%
        promoter_conn_matrix[,distal_safe_sites, drop=FALSE] %*%
        scaled_site_weights
    scaled_site_weights@x[scaled_site_weights@x > 1] <- 1

    # Multiply distal accessibility by site weights
    distal_activity_scores <- scaled_site_weights %*% distal_activity_scores

    dimnames(distal_activity_scores) <- dimnames(promoter_access_mat_in_cicero_map)


    distal_activity_scores <-
        distal_activity_scores[row.names(promoter_access_mat_in_cicero_map),, drop=FALSE]

    # Sum distal and promoter scores
    promoter_activity_scores <- distal_activity_scores +
        promoter_access_mat_in_cicero_map

    # Make and populate final matrix
    promoter_gene_mat <-
        Matrix::sparseMatrix(j=as.numeric(factor(promoter_peak_table$peak)),
                             i=as.numeric(factor(promoter_peak_table$gene)),
                             x=1)
    colnames(promoter_gene_mat) = levels(factor(promoter_peak_table$peak))
    row.names(promoter_gene_mat) = levels(factor(promoter_peak_table$gene))
    promoter_gene_mat <- promoter_gene_mat[,row.names(promoter_activity_scores)]
    gene_activity_scores <- promoter_gene_mat %*% promoter_activity_scores

    return(gene_activity_scores)
}




###############################################
######## TO EDIT FOR MY USE CASE
################################################

coaccess <- read.table(coaccess)

# foreach parameters
message(" ... Initializing per cluster cicero run - GENE ACTIVITY")
gascores <- lapply(clusts, function(x){

    # get gene activity scores --------------------------------------------------------------------
    message("--- estimating gene activity scores for cluster ",x)
    ids <- rownames(meta2[meta2$Cluster==x,])
    index.keep <- colnames(exprs(cds)) %in% ids
    s.cds <- cds[,index.keep]
    
    keep_a <- rownames(coaccess[coaccess$group==x,])
    a.sub <- coaccess[keep_a,]
    a.sub <- a.sub[,c(1,2,3)]
    colnames(a.sub) <- c("Peak1","Peak2","coaccess")
    
    # only consider sites accessible in at least 1% of cells in cluster
    s.cds <- s.cds[Matrix::rowSums(exprs(s.cds))>0,]
    s.cds <- s.cds[,Matrix::colSums(exprs(s.cds))>0]
    print(head(exprs(s.cds)[,1:5]))
    message(" - number of sites for cluster ", x, " = ", nrow(s.cds))
    
    # get UMAP coordinates
    umap_coords <- pData(s.cds)[,c("umap1","umap2")]

    message("# UMAP coords = ", nrow(umap_coords), " | # cells = ", ncol(s.cds))
    cicero_cds <- make_cicero_cds(s.cds, reduced_coordinates=umap_coords, k=30)


    
    # estimate gene activity
    message(" ... Estimating gene activity scores")
    pos <- subset(genes, strand == "+")
    pos <- pos[order(pos$start),] 
    pos <- pos[!duplicated(pos$transcript),]
    pos$end <- pos$start
    pos$start <- pos$start - 1000
    neg <- subset(genes, strand == "-")
    neg <- neg[order(neg$start, decreasing = T),] 
    neg <- neg[!duplicated(neg$transcript),] 
    neg$start <- neg$end
    neg$end <- neg$end + 1000
    
    # merge
    gene_ann2 <- rbind(pos, neg)
    gene_ann2 <- gene_ann2[,c(1:3, 5)]
    gene_ann2 <- gene_ann2[order(gene_ann2$start, decreasing=F),]
    colnames(gene_ann2)[4] <- "gene"
    
    # annotate genes
    message("     - annotate genes by peaks ...")
    s.cds <- annotate_cds_by_site(s.cds, gene_ann2, all=F,header=T)
    
    # estimate un-normalized activity
    message("     - build gene activity matrix ... ")
    saveRDS(s.cds, file=paste("s.cds_cluster",x,".rds",sep=""))
    saveRDS(a.sub, file=paste("a.sub_cluster",x,".rds",sep=""))

    dist_thresh=250000
    coaccess_cutoff=0.25

    accessibility_mat <- exprs(s.cds)
    site_weights <- Matrix::rowMeans(accessibility_mat)/Matrix::rowMeans(accessibility_mat)
    site_weights[names(site_weights)] <- 1

    input_cds <- s.cds

    cicero_cons_info <- a.sub

    unnorm_ga <- build_composite_gene_activity_matrix(input_cds, site_weights, cicero_cons_info,dist_thresh = dist_thresh,coaccess_cutoff = coaccess_cutoff)

    write.table(unnorm_ga, file=paste(out,"Cluster",x,".","unnormalized.cicero.geneActivity.txt",sep=""),
                sep="\t",quote=F, col.names=T, row.names=T)
    unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, !Matrix::colSums(unnorm_ga) == 0]
    


    # gene activity per cluster
    message("      gene activity per cluster... ")
    num_genes <- pData(s.cds)$num_genes_expressed
    names(num_genes) <- row.names(pData(s.cds))


    
    
    # normalize
    message("      normalize gene activities ... ")
    cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)
    geneact <- as.data.frame(summary(cicero_gene_activities))
    geneact$i <- rownames(cicero_gene_activities)[geneact$i]
    geneact$j <- colnames(cicero_gene_activities)[geneact$j]
    
    # output
    write.table(geneact, file=paste(out,"Cluster",x,".","cicero.geneActivity.txt",sep=""),
                sep="\t",quote=F, col.names=F, row.names=F)
    
    # return
    return(unnorm_ga)
})
