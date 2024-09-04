#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=1-00:15:00
#SBATCH --job-name="004_FRIP_QC.sh"
#SBATCH --output=std/004_FRIP_QC%j.stdout
#SBATCH --output=std/004_FRIP_QC%j.stderr
#SBATCH -p seymourlab
#SBATCH --array=1-21


conda activate Deeptools
module load bedtools
module load samtools

SEQLIST=/bigdata/seymourlab/idiaz026/Data/Floral_scATAC/2021-07-03_20WGS_samplelist.txt
BAMS=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-05-05_Downsample_ACR_analysis/DownsampledBams/Clean_Bams
PEAKS=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_by_geno/Downsampled_genrich
COLLATE_PEAKS=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-05-05_Downsample_ACR_analysis/DownsampledBams/Peaks/All_ATAC.peaks.bed

PLOTS=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-05-05_Downsample_ACR_analysis/FRIP

NAME=$(head -n "$SLURM_ARRAY_TASK_ID" "$SEQLIST" | tail -n1 | cut -f3)
ID=$(head -n "$SLURM_ARRAY_TASK_ID" "$SEQLIST" | tail -n1 | cut -f1)

TEMP_DIR=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-05-05_Downsample_ACR_analysis/Temp
TEMP_DIR="$TEMP_DIR"/"$NAME"_"$RANDOM"
mkdir -pv "$TEMP_DIR"

cd $TEMP_DIR

samtools index "$BAMS"/"$ID"_rmdup_mq10.bam

#All_samples at once

#plotEnrichment --bamfiles $(ls "$BAMS"/*_rmdup_mq10.bam | tr '\n' ' ') --BED $COLLATE_PEAKS --plotFile $PLOTS/All_Genotypes.collate.orig.png --labels $(ls "$BAMS"/*mq10.bam |cut -c 91-| tr '\n' ' ') --outRawCounts "$PLOTS"/Collate_peaks_allFrip.txt


plotEnrichment --bamfiles "$BAMS"/"$ID"_rmdup_mq10.bam --BED $PEAKS/"$ID"*CLEAN.bed --plotFile $PLOTS/$ID.orig.png -p max --labels $ID --outRawCounts "$PLOTS"/$ID.orig.frip.txt
