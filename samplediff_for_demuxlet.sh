#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=32G
#SBATCH --output=std/samplediff_for_demuxlet.stderr
#SBATCH --error=std/samplediff_for_demuxlet.stdout
#SBATCH --job-name="samplediff_for_demuxlet.sh"
#SBATCH --time=5-00:00:00
#SBATCH -p seymourlab


source activate /rhome/idiaz026/.conda/envs/Floral_scATAC

VCF=/bigdata/seymourlab/idiaz026/Results/Floral_scATAC/vcf/Merged_filtering/2021-09-25_99_miss_filt3_reheader.vcf.gz

cd /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2022-08-25_process_scATAC/Demuxlet/sample_diffs

plink2 --vcf $VCF --aec --out VCF_99miss.diff --sample-diff pairwise id-delim=, ids=Tango Valentine Washington_navel Algerian_clementine USDA_88-2 Koster_mandarin Frost_owen_satsuma Oroblanco_grpft Orlando Bouquettier_de_Nice Citrus_junos Seedless_kishu Interdonato_lemon Frost_lisbon S.emperor_ponkan Kao_panne_pummelo Swingle SCFS_citron Microcitrus_australiasica C_macrophylla

echo "Based on most similar sample (manually determined), select sdiff files to extract informative SNPs"

cat VCF_99miss.diff.Tango.S.emperor_ponkan.sdiff VCF_99miss.diff.Valentine.Oroblanco_grpft.sdiff VCF_99miss.diff.Washington_navel.Orlando.sdiff VCF_99miss.diff.Algerian_clementine.Koster_mandarin.sdiff VCF_99miss.diff.USDA_88-2.Orlando.sdiff VCF_99miss.diff.Koster_mandarin.S.emperor_ponkan.sdiff VCF_99miss.diff.Frost_owen_satsuma.Seedless_kishu.sdiff VCF_99miss.diff.Washington_navel.Bouquettier_de_Nice.sdiff VCF_99miss.diff.Citrus_junos.Seedless_kishu.sdiff VCF_99miss.diff.Interdonato_lemon.Frost_lisbon.sdiff VCF_99miss.diff.Oroblanco_grpft.Kao_panne_pummelo.sdiff VCF_99miss.diff.Kao_panne_pummelo.Swingle.sdiff VCF_99miss.diff.Kao_panne_pummelo.Microcitrus_australiasica.sdiff VCF_99miss.diff.Kao_panne_pummelo.C_macrophylla.sdiff | cut -f1,2 | sort -k1,1 -k2,2n | uniq | tail -n +2  > most_informative_SNPs.txt



