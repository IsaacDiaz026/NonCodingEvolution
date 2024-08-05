#!/usr/bin/env Rscript

library(devtools)

library(Socrates)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

bed <- args[1]

ann <- "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/WN_HAPA_genes_modCTG.gff.gz"

chr <- "/bigdata/seymourlab/idiaz026/REFERENCES/PWN_HAPA/WN_HAPA_modContig_genomeFile.txt"

obj <- loadBEDandGenomeData(bed, ann, chr, is.fragment =T)


print("Peak calling")
obj <- callACRs(obj, genomesize=3.0e8,
                shift= -50,
                extsize=100,
                fdr=0.05,
                output=paste(args[2],"_postClust",sep=""),
                tempdir=args[2],
                verbose=T)




