#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=700G
#SBATCH --job-name="000_cellranger_atac_count.sh"
#SBATCH --output=std/000_cellranger_atac_count.sh.out
#SBATCH --time=12-00:15:00
#SBATCH -p highmem

cd /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-12_cellrangerPWN_redo

TMPDIR=/rhome/idiaz026/.tmp
GFF=/bigdata/seymourlab/eavilade/NLR/Genomes/WN_HapA.fa.gb_1.1c_2.3_s1_0_16_1.gff
FASTQS=/bigdata/seymourlab/idiaz026/Data/Floral_scATAC/2023-02-20_raw_scATAC_fastq


#/bigdata/seymourlab/idiaz026/Apps/cellranger-atac-2.1.0/cellranger-atac mkref --config=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-12_cellrangerPWN_redo/PWN.config

/bigdata/seymourlab/idiaz026/Apps/cellranger-atac-2.1.0/cellranger-atac count --id=PWN_HAPA_wPlastid_fix --reference=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-12_cellrangerPWN_redo/PWNplusPlastid --fastqs /bigdata/seymourlab/idiaz026/Data/Floral_scATAC/2023-02-20_raw_scATAC_fastq
