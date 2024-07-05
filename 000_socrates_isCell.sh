#!/bin/bash -l

#SBATCH --ntasks=16
#SBATCH --mem=900G
#SBATCH --job-name="000_socrates_isCell.sh"
#SBATCH --output=std/000_socrates_isCell.sh.out
#SBATCH --time=12-00:15:00
#SBATCH -p highmem


conda activate Socrates_fix

ALL_ATAC=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-12_cellrangerPWN_redo/PWN_HAPA_wPlastid_fix/outs/possorted_bam.bam
TMPDIR=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Temp

INPUT=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-12_cellrangerPWN_redo/PWN_HAPA_wPlastid_fix/outs

cd /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates

#sinto fragments -b "$ALL_ATAC" -f "$INPUT"/All_ATAC_plastidfragments.bed -p 16 --use_chrom "(?i)^NC-"

#gunzip "$INPUT"/All_ATAC_plastidfragments.bed
#gunzip "$INPUT"/fragments.tsv

#cat "$INPUT"/fragments.tsv "$INPUT"/All_ATAC_plastidfragments.bed > "$INPUT"/All_ATAC_plusplastidfragments.bed

#gzip "$INPUT"/All_ATAC_plusplastidfragments.bed

#gzip /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Tn5_frag_conv2bed.bed

#export R_MAX_VSIZE=900Gb
Rscript 000_socrates_isCell.R

gzip All_ATAC.sparse

cut -f1,18 2023-08-24_raw_cell_data_isCELL.txt | tail -n +2 > barcodes_to_keep.txt

sinto filterbarcodes -b $ALL_ATAC -c barcodes_to_keep.txt -p 16 --outdir /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Demuxlet

mv /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Demuxlet/A.bam /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Demuxlet/scATAC_filtered.bam

