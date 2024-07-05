#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=80G
#SBATCH --time=2-00:15:00
#SBATCH --job-name="000_makeGenoFragFiles.sh"
#SBATCH --output=std/000_makeGenoFragFiles.sh.sh%j.stdout
#SBATCH --output=std/000_makeGenoFragFiles.sh%j.stderr
#SBATCH -p intel
#SBATCH --array=1-21

conda activate Socrates
module load bedtools
SEQLIST=/bigdata/seymourlab/idiaz026/Data/Floral_scATAC/2021-07-03_20WGS_samplelist.txt
NAME=$(head -n "$SLURM_ARRAY_TASK_ID" "$SEQLIST" | tail -n1 | cut -f3)
ID=$(head -n "$SLURM_ARRAY_TASK_ID" "$SEQLIST" | tail -n1 | cut -f1)

BAM=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/SintoGenotypeBams/Collate_BAM/All_ATAC.clean.bam
FRAGMENTS=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/SintoGenotypeBams/scATACFragments_files
PEAKS=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_by_geno

TMPDIR=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Temp

GENOME_FILE=/bigdata/seymourlab/idiaz026/REFERENCES/PWN_HAPA/WN_HAPA_modContig_genomeFile.txt

cd /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/

#sinto filterbarcodes -b $BAM -c "$FRAGMENTS"/"$ID".postClust.barcodes.txt --outdir "$FRAGMENTS"

#samtools index "$FRAGMENTS"/"$ID".bam
#sinto fragments -b "$FRAGMENTS"/"$ID".bam -f $FRAGMENTS/"$ID"_fragments.bed --use_chrom "(?i)^Scaffold-" --collapse_within 

#sort -k 1,1 -k2,2n $FRAGMENTS/"$ID"_fragments.bed > $FRAGMENTS/"$ID"_fragments.sort.bed

#bgzip -@ 8 $FRAGMENTS/"$ID"_fragments.sort.bed
#tabix -p bed $FRAGMENTS/"$ID"_fragments.sort.bed.gz

Rscript 003c_peak_calling.R  $FRAGMENTS/"$ID"_fragments.sort.bed.gz "$PEAKS"/"$ID"









