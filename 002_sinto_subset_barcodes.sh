#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=600G
#SBATCH --time=12-00:15:00
#SBATCH --job-name="002_sinto_subset_barcodes.sh"
#SBATCH --output=std/002_sinto_subset_barcodes.sh.sh.out
#SBATCH -p highmem

conda activate Socrates

ALL_ATAC=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-12_cellrangerPWN_redo/PWN_HAPA/outs/possorted_bam.bam

cd /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/SintoGenotypeBams

sinto filterbarcodes -b $ALL_ATAC -c /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/SintoGenotypeBams/2023-09-06_all_snps_barcodesClean.txt --outdir /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/SintoGenotypeBams

cd /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/SintoGenotypeBams/

BAMS=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/SintoGenotypeBams/Collate_BAM
TMPDIR=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Temp

samtools merge $BAMS/All_ATAC.bam *.bam

samtools index $BAMS/All_ATAC.bam

java -jar /opt/linux/centos/7.x/x86_64/pkgs/picard/2.18.3/lib/picard.jar MarkDuplicates INPUT= $BAMS/All_ATAC.bam O= $BAMS/All_ATAC.rmdup.bam M= $BAMS/All_ATAC_rmdup_metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=false VALIDATION_STRINGENCY=LENIENT BARCODE_TAG=CB

samtools index $BAMS/All_ATAC.rmdup.bam

samtools view -f 3 -bhq 10 $BAMS/All_ATAC.rmdup.bam > "$BAMS"/All_ATAC.clean.bam

samtools index All_ATAC.clean.bam

sinto fragments -b "$BAMS"/All_ATAC.clean.bam -f "$BAMS"/All_ATAC_fragments.bed --use_chrom "(?i)^Scaffold-" --collapse_within

gzip "$BAMS"/All_ATAC_fragments.bed




#sinto filterbarcodes -b $ALL_ATAC -c inform_snps_barcodesClean.txt --outdir /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2022-08-25_process_scATAC/Demuxlet/2023-04-01_sintoExtractBarcodes






