#2023-03-24 plotting Demuxlet output


setwd("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Demuxlet/")
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(dplyr)

Floral_scATAC_fullBam <- read.delim("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Demuxlet/Floral_scATAC_fullBam.best")
Floral_scATAC_fullBam$full_geno <- Floral_scATAC_fullBam$SNG.BEST.GUESS


Floral_scATAC <- Floral_scATAC_fullBam


#assign type based on difference between likelihoods 

Floral_scATAC$LRT.SNG <- Floral_scATAC$SNG.BEST.LLK - Floral_scATAC$SNG.NEXT.LLK

Floral_scATAC$type <- ifelse(Floral_scATAC$LRT.SNG >1 , "singlet","doublet")

table(Floral_scATAC$type)

singlets <- Floral_scATAC %>% filter(type=="singlet")

ggplot(singlets, aes(x = SNG.BEST.GUESS, y =LRT.SNG)) +  stat_boxplot(geom = "errorbar", # Error bars
              width = 0.25) +    # Bars width
 geom_boxplot() + coord_flip() + theme_classic() + ylim(0,50)


ggplot(singlets, aes(x = SNG.BEST.GUESS, y =BEST.POSTERIOR)) + 
 stat_boxplot(geom = "errorbar", # Error bars
              width = 0.25) +    # Bars width
 geom_boxplot() + coord_flip() + theme_classic()

#######################
Floral_scATAC$LRT.SNG <- Floral_scATAC$SNG.BEST.LLK - Floral_scATAC$SNG.NEXT.LLK

Floral_scATAC$type <- ifelse(Floral_scATAC$LRT.SNG >1 , "singlet","doublet")

table(Floral_scATAC$type)

singlets <- Floral_scATAC %>% filter(type=="singlet")

ggplot(singlets, aes(x = SNG.BEST.GUESS, y =LRT.SNG)) +  stat_boxplot(geom = "errorbar", # Error bars
                                                                      width = 0.25) +    # Bars width
 geom_boxplot() + coord_flip() + theme_classic() + ylim(0,50)


#####################
############################


best_barcodes <- singlets[,c(2,21)]


#import barcodes that correspond to cells called by "is_cell"
barcodes_to_keep <- read.delim("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/barcodes_to_keep.txt", header=FALSE)

true_barcodes <- best_barcodes[(best_barcodes$BARCODE %in% barcodes_to_keep$V1),]


table(true_barcodes$full_geno)



Summary <- read.table("/bigdata/seymourlab/idiaz026/Data/Floral_scATAC/scATAC_stats_open_chrom.txt", header = T, sep = "\t")

add_geno <- Summary[,c(1,2)]


true_barcodes <- merge(true_barcodes, add_geno, by.x = "full_geno",by.y="Sample_ID")

b <- table(true_barcodes$genotype)

write.table(b,"2023-09-06_NumNuc_pergeno.txt",sep="\t",quote=F, row.names=F, col.names=F)



true_barcodes <- true_barcodes[,c(2,3)]

write.table(true_barcodes, "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/SintoGenotypeBams/2023-09-06_all_snps_barcodesClean.txt",sep = "\t",col.names = F, row.names = F, quote=F)

minus_inter_cmac <- true_barcodes %>% filter(genotype != "Interdonato_lemon")

write.table(minus_inter_cmac, "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/SintoGenotypeBams/2023-09-06_all_snps_barcodesClean_minusInter.txt",sep = "\t",col.names = F, row.names = F, quote=F)



