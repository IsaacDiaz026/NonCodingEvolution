#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=600G
#SBATCH --time=9-00:15:00
#SBATCH --job-name="Demuxlet_scATAC.sh"
#SBATCH --output=std/Demuxlet_scATAC.stdout
#SBATCH --output=std/Demuxlet_scATAC.stderr
#SBATCH -p highmem


conda activate Demuxlet
module load samtools
module load bcftools

VCF=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-07-04_WGS_variant_callingPWN/vcf/2023-07-04_filter3.vcf.gz

BAM=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-12_cellrangerPWN_redo/PWN_HAPA/outs/possorted_bam.bam

OUTPUT=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Demuxlet


cd /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates


#bcftools annotate --rename-chrs rename_ctgs.txt $VCF > 2023-08-20_scATAC_wgs_renamed.vcf

#bgzip 2023-08-20_scATAC_wgs_renamed.vcf

#tabix -p vcf 2023-08-20_scATAC_wgs_renamed.vcf.gz

VCF=2023-08-20_scATAC_wgs_renamed.vcf.gz
popscle demuxlet --sam $BAM --vcf $VCF --field GT --out $OUTPUT/Floral_scATAC_fullBam
