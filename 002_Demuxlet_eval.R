#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(dplyr)
library(argparse)



args <- commandArgs(trailingOnly=TRUE)
#load demuxlet output
Floral_scATAC <- fread(args[1])

#import barcodes that correspond to cells called by Socrates "is_cell"
barcodes_to_keep <- fread(args[2])

#sample names with genotype info
add_geno <- fread(args[3])



Floral_scATAC$full_geno <- Floral_scATAC$SNG.BEST.GUESS

#assign type based on difference between likelihoods 
Floral_scATAC$LRT.SNG <- Floral_scATAC$SNG.BEST.LLK - Floral_scATAC$SNG.NEXT.LLK # nolint

Floral_scATAC$type <- ifelse(Floral_scATAC$LRT.SNG > 1 , "singlet","doublet")

table(Floral_scATAC$type)

singlets <- Floral_scATAC %>% filter(type=="singlet")

pdf("DEMUX/demuxlet_report.pdf")
ggplot(singlets, aes(x = SNG.BEST.GUESS, y =LRT.SNG)) + 
    stat_boxplot(geom = "errorbar",width = 0.25) +
    geom_boxplot() + coord_flip() + theme_classic() +
    ylim(0, 50)


ggplot(singlets, aes(x = SNG.BEST.GUESS, y =BEST.POSTERIOR)) +
 stat_boxplot(geom = "errorbar", # Error bars
              width = 0.25) +    # Bars width
 geom_boxplot() + coord_flip() + theme_classic()

#######################
Floral_scATAC$LRT.SNG <- Floral_scATAC$SNG.BEST.LLK - Floral_scATAC$SNG.NEXT.LLK

Floral_scATAC$type <- ifelse(Floral_scATAC$LRT.SNG >1 , "singlet","doublet")

table(Floral_scATAC$type)

singlets <- Floral_scATAC %>% filter(type=="singlet")

ggplot(singlets, aes(x = SNG.BEST.GUESS, y =LRT.SNG)) +
    stat_boxplot(geom = "errorbar",width = 0.25) +
    geom_boxplot() + coord_flip() + theme_classic() + ylim(0, 50)

dev.off()


best_barcodes <- singlets[,c(2,21)]

#clean barcodes using info from socrates is cell
true_barcodes <- best_barcodes[(best_barcodes$BARCODE %in% barcodes_to_keep$V1),]


table(true_barcodes$full_geno)


true_barcodes <- merge(true_barcodes, add_geno, by.x = "full_geno", by.y = "Sample_ID")

b <- table(true_barcodes$genotype)

write.table(b, "DEMUX/NumNuc_pergeno.txt",sep="\t",
            quote = FALSE, row.names = FALSE, 
            col.names = FALSE)


true_barcodes <- true_barcodes[,c(2,3)]

write.table(true_barcodes, "DEMUX/Nuclei_barcodesClean.txt", sep = "\t",
            col.names = FALSE, row.names = FALSE, quote = FALSE)

#remove Interdonato
minus_interdonato <- true_barcodes %>% filter(genotype != "Interdonato_lemon")

write.table(minus_inter, "DEMUX/barcodesClean_minusInter.txt",
            sep = "\t", col.names = FALSE,
            row.names = FALSE, quote = FALSE)