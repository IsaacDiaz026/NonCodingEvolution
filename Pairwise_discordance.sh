#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=32G
#SBATCH --output=std/005_VCF_pca.stderr
#SBATCH --error=std/005_VCF_pca.stdout
#SBATCH --job-name="005_VCF_pca.sh"
#SBATCH --time=5-00:00:00
#SBATCH -p seymourlab


source activate /rhome/idiaz026/.conda/envs/Floral_scATAC
module load bcftools
module load plink



SAMPLELIST=/bigdata/seymourlab/idiaz026/Data/Floral_scATAC/2021-09-07_20WGS_samplelist.txt

VCF=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-07-04_WGS_variant_callingPWN/vcf/2023-07-04_filter3.vcf.gz


RESULTS=/bigdata/seymourlab/idiaz026/NonCodingEvolution/Plink

cd $RESULTS




ALL_ACRs=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/ACRs_classified.sort.bed

cat $ALL_ACRs | cut -f1,2,3 | sed 's/Scaffold-/Scaffold_/g' > ALL_acrs.bed

ALL_ACRs=ALL_acrs.bed

bcftools view -O v -R $ALL_ACRs $VCF > VCF_75.only_ACRs.vcf


plink2 --vcf VCF_75.only_ACRs.vcf --aec --out VCF_75_miss.diff --sample-diff id-delim=, ids=IAD_01_S367 IAD_02_S368 IAD_03_S369 IAD_04_S370 IAD_05_S371 IAD_06_S372 IAD_07_S373 IAD_08_S374 IAD_09_S375 IAD_10_S376 IAD_11_S377 IAD_12_S378 IAD_14_S380 IAD_15_S381 IAD_16_S382 IAD_17_S383 IAD_18_S384 IAD_19_S385 IAD_20_S386




