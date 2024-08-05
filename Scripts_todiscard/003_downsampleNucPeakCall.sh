#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --mem=80G
#SBATCH --output=std/000_downsampleBarcodesPeakCall%j.stdout
#SBATCH --job-name="000_downsampleBarcodesPeakCall%j.sh"
#SBATCH -p intel
#SBATCH --array=1-21



source activate /rhome/idiaz026/.conda/envs/Socrates
module load samtools


ALL_ATAC=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/SintoGenotypeBams/Collate_BAM/All_ATAC.clean.bam

BARCODES=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/SintoGenotypeBams/2023-09-06_all_snps_barcodesClean_minusInter.txt

SEQLIST=/bigdata/seymourlab/idiaz026/Data/Floral_scATAC/2021-07-03_20WGS_samplelist.txt

GENO_BARCODES=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-05-05_Downsample_ACR_analysis/GenotypeExtractedBarcodes

BAMS=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-05-05_Downsample_ACR_analysis/DownsampledBams

CLEAN_BAMS=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-05-05_Downsample_ACR_analysis/DownsampledBams/Clean_Bams

DOWNSAMPLE_BARCODES=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-05-05_Downsample_ACR_analysis/GenotypeExtractedBarcodes/DownsampledBarcodes

TEMP_DIR=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-05-05_Downsample_ACR_analysis/Temp

RM_DUP=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-05-05_Downsample_ACR_analysis/DownsampledBams/rmdup

INPUT=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-07-04_WGS_variant_callingPWN/Alignments/Merged_bams/Correct_RG/Mark_dup

RESULTS=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-05-05_Downsample_ACR_analysis/DownsampledBams/Peaks

PCR_DUP=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-05-05_Downsample_ACR_analysis/DownsampledBams/PCR_DUP

TEMP_DIR="$TEMP_DIR"/"$NAME"_"$RANDOM"
mkdir -pv "$TEMP_DIR"

NAME=$(head -n "$SLURM_ARRAY_TASK_ID" "$SEQLIST" | tail -n1 | cut -f3)
ID=$(head -n "$SLURM_ARRAY_TASK_ID" "$SEQLIST" | tail -n1 | cut -f1)


cat $BARCODES | grep "$ID" > "$GENO_BARCODES"/"$ID".barcodes.txt

num_barcodes=$(cat "$GENO_BARCODES"/"$ID".barcodes.txt | wc -l)


echo "If number of barcodes is greater than 552, downsample 552 barcodes randomly. 552 is median number of barcodes across 20 genotypes"

if [[ $num_barcodes -gt 552 ]]
then
    cat $GENO_BARCODES/"$ID".barcodes.txt | shuf -n 552 > "$DOWNSAMPLE_BARCODES"/"$ID"_ds_barcodes.txt
else
    cat $GENO_BARCODES/"$ID".barcodes.txt > "$DOWNSAMPLE_BARCODES"/"$ID"_ds_barcodes.txt
fi


sinto filterbarcodes -b $ALL_ATAC -c "$DOWNSAMPLE_BARCODES"/"$ID"_ds_barcodes.txt -p 8 --outdir /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-05-05_Downsample_ACR_analysis/DownsampledBams

samtools index $BAMS/"$ID".bam

echo "Removing duplicates from bam"
#java -jar /opt/linux/centos/7.x/x86_64/pkgs/picard/2.18.3/lib/picard.jar MarkDuplicates INPUT= $BAMS/$ID.bam O= $RM_DUP/$ID.rmdup.bam M= $RM_DUP/metrics/"$ID"rmdup_metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT BARCODE_TAG=CB


#samtools index $RM_DUP/$ID.rmdup.bam

#samtools view -f 3 -bhq 10 $RM_DUP/$ID.rmdup.bam > $CLEAN_BAMS/"$ID"_rmdup_mq10.bam

samtools sort -n $BAMS/"$ID".bam -o $TEMP_DIR/$ID.namesort.bam

samtools sort -n $INPUT/"$ID".marked_duplicates.bam -o $TEMP_DIR/$ID.namesort.input.bam


~/Apps/Genrich-master/Genrich -t $TEMP_DIR/$ID.namesort.bam -R $PCR_DUP/$ID.pcrdups.txt -p 0.05 -r -y -j -o $RESULTS/$ID.05.peaks.txt -f $RESULTS/$ID.05.signal.log -c $TEMP_DIR/$ID.namesort.input.bam

rm -r $TEMP_DIR












