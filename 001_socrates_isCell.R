#!/usr/bin/env Rscript

library(devtools)
library(Socrates)
library(dplyr)
library(doSNOW)
library(argparse)

args <- commandArgs(trailingOnly=TRUE)


bed <- args[1]
ann <- args[2]
chr <- args[3]

macsPath=args[4]
orgToFilt=args[5]


obj <- loadBEDandGenomeData(bed, ann, chr, is.fragment =T)

print("print_obj_size")

write.table(obj$bed, "Tn5_frag_conv2bed.bed", sep = "\t", col.names=FALSE, row.names=FALSE,quote=FALSE)

obj <- countRemoveOrganelle(obj, org_scaffolds=c("NC-037463","NC-034671"), remove_reads=FALSE)           

obj <- callACRs(obj, genomesize=3.0e8,
                shift= -75,
                extsize=100,
                fdr=0.05,
                output="bulk_peaks",
                tempdir=paste(macsPath),
                verbose=T)

print("build meta data")

obj <- buildMetaData(obj, tss.window=2000, verbose=TRUE, 
    organelle_scaffolds =c("NC-037463","NC-034671"))
#system(paste0("cat /proc/",Sys.getpid(),"/status | grep VmSize"))
print("head obj after buildMetaData")
head(obj$meta)

obj <- findCells(obj,
                 doplot=T,
                 min.cells=1000,
                 max.cells=10000,
                 min.tn5=1000,
		 filt.org=T,
                 org.filter.thresh=0.8,
                 tss.min.freq = 0.1,
                 tss.z.thresh=3,
                 filt.tss=TRUE,
                 filt.frip=TRUE,
                 frip.min.freq=0.1,
                 prefix = "scATAC_qualC")

print("generate matrix")
obj <- generateMatrix(obj, filtered=F, 
    peaks=F,windows=500, verbose=T)


write.table(obj$counts, 
    "All_ATAC.sparse", 
    row.names= F, col.names = F, quote = F , sep = "\t")

print("head obj$meta after gen matrix")
print(head(obj$meta))
print("dim obj$meta after gen matrix")
print(dim(obj$meta))


print("dim of obj$counts after gen matrix")
print(dim(obj$counts))
print("head obj$counts after gen matrix")
print(head(obj$counts))

write.table(obj$meta, "cell_meta_beforeisCell.txt", col.names=T, row.names=T, sep = "\t", quote=FALSE)

print("running is cell")
obj <- isCell(obj,verbose = T,min.FRiP = 0.1)

print("running is cell complete")

print("dim obj$meta after isCell")
print(dim(obj$meta))

print("head obj$meta")
head(obj$meta)

print("make df")
raw_cell_data <- as.data.frame(obj$meta)
print(head(raw_cell_data))
print(summary(raw_cell_data$is_cell))
print(summary(raw_cell_data$background))

write.table(raw_cell_data,"raw_cell_data_beforeFilt.txt",col.names=T, 
    row.names=F, sep = "\t",quote=FALSE)


print("filt df")
raw_cell_data <- raw_cell_data %>% filter(is_cell == 1)
print("remake_df")
raw_cell_data <- as.data.frame(raw_cell_data)
print(class(raw_cell_data))
print("add column")
raw_cell_data$group <- "A"


write.table(raw_cell_data,"raw_cell_data_isCELL.txt",
    col.names=T, row.names=F, sep = "\t",quote=FALSE)


soc.obj <- convertSparseData(obj, verbose=T)


saveRDS(obj, file="QC_object.rds")
