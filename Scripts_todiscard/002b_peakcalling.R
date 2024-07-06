#!/usr/bin/env Rscript

library(devtools)

library(Socrates)
library(dplyr)
library(doSNOW)
library(ggplot2)
library(stringr)

#########################
##########EDIT BED FOR CLEAN FRAG FILE####
bed <- "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/SintoGenotypeBams/Collate_BAM/All_ATAC_fragments.bed.gz"
##############

ann <- "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/WN_HAPA_genes_modCTG.gff.gz"

chr <- "/bigdata/seymourlab/idiaz026/REFERENCES/PWN_HAPA/WN_HAPA_modContig_genomeFile.txt"

obj <- loadBEDandGenomeData(bed, ann, chr, is.fragment =T)

write.table(obj$bed, "Clean_Tn5_frag_conv2bed.bed", sep = "\t", col.names=FALSE, row.names=FALSE,quote=FALSE)

print("Peak calling")
obj <- callACRs(obj, genomesize=3.0e8,
                shift= -50,
                extsize=100,
                fdr=0.05,
                output="bulk_peaks_CLEAN",
                tempdir="./macs2_temp",
                verbose=T)

obj <- buildMetaData(obj, tss.window=2000, verbose=TRUE)


print("Generate matrix using only peaks")
obj <- generateMatrix(obj, filtered=F, windows=1000, peaks=T, verbose=T)

print("dim obj $ counts")
print(dim(obj$counts))

write.table(obj$counts,"/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/All_ATAC_peaks.CLEAN.sparse", row.names= F, col.names = F, quote = F , sep = "\t")

soc.obj <- convertSparseData(obj, verbose=T)



peaks <- obj$counts 

geno <- read.table("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/SintoGenotypeBams/2023-08-27_OLD_isCell_all_snps_barcodesClean_minusInterCmac.txt")

colnames(geno) <- c("CellID","Genotype")

colnames(peaks) <- c("Region", "CellID","Open")

peaks$chr <- str_split_i(peaks$Region, "_",1)

peaks$start <- str_split_i(peaks$Region, "_",2)

peaks$end <- str_split_i(peaks$Region, "_",3)

print("Attaching genotype information")

peaks <- merge(peaks, geno, by=c("CellID"))

genotypes <- unique(peaks$Genotype)

dfs <- list()
for (i in 1:length(genotypes)){
    dfs[[i]] <- peaks %>% filter(Genotype == genotypes[i])
    dfs[[i]] <- 
    write.table(dfs[[i]][,c(4,5,6)], file=paste("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_by_geno/",genotypes[i],"_ACRs.bed",sep=""),sep="\t",quote=F, row.names=F, col.names=F)
    
}















