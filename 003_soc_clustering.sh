#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mem=600G
#SBATCH --job-name="003_soc_clustering_strict.sh"
#SBATCH --output=std/003_soc_clustering_strict.sh.out
#SBATCH -p highmem

conda activate Socrates
module load bedtools

#Rscript /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/003_soc_clustering.R

GENES=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Counts_Gene_Activity/Genes_upAndDownstream.bed

BARCO#DES=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/2023-09-09_cellbarcodes_clean.txt

TN5=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Clean_Tn5_frag_conv2bed.bed.gz

TMPDIR=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Temp
cd /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Counts_Gene_Activity/

BAM=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/SintoGenotypeBams/Collate_BAM/All_ATAC.clean.bam

OUT=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Counts_Gene_Activity

Clusters=(Cluster1 Cluster2 Cluster3 Cluster4)


#for i in ${!Clusters[@]}; do cut -f1,12 ${Clusters[i]}_Rep1.metadata.txt > ${Clusters[i]}_Rep1.barcodes.txt
#cut -f1,12 ${Clusters[i]}_Rep2.metadata.txt > ${Clusters[i]}_Rep2.barcodes.txt
#cut -f1,12 ${Clusters[i]}_Rep3.metadata.txt > ${Clusters[i]}_Rep3.barcodes.txt  
#sinto filterbarcodes -b $BAM -c ${Clusters[i]}_Rep1.barcodes.txt --outdir "$OUT"
#sinto filterbarcodes -b $BAM -c ${Clusters[i]}_Rep2.barcodes.txt --outdir "$OUT"
#sinto filterbarcodes -b $BAM -c ${Clusters[i]}_Rep3.barcodes.txt --outdir "$OUT"

#bedtools coverage -a "$GENES" -b "${Clusters[i]}"_rep1.bam > /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Counts_Gene_Activity/Adjusted/Gene_counts_"${Clusters[i]}"_Rep1.txt
#bedtools coverage -a "$GENES" -b "${Clusters[i]}"_rep2.bam > /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Counts_Gene_Activity/Adjusted/Gene_counts_"${Clusters[i]}"_Rep2.txt
#bedtools coverage -a "$GENES" -b "${Clusters[i]}"_rep3.bam > /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Counts_Gene_Activity/Adjusted/Gene_counts_"${Clusters[i]}"_Rep3.txt
#done
#Clusters=(Cluster_1 Cluster_2 Cluster_3 Cluster_4)

#for i in ${!Clusters[@]}; do bedtools coverage -a "$GENES" -b "${Clusters[i]}"_rep1.bam > /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Counts_Gene_Activity/Adjusted/Gene_counts_"${Clusters[i]}"_Rep1.txt
#bedtools coverage -a "$GENES" -b "${Clusters[i]}"_rep2.bam > /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Counts_Gene_Activity/Adjusted/Gene_counts_"${Clusters[i]}"_Rep2.txt
#bedtools coverage -a "$GENES" -b "${Clusters[i]}"_rep3.bam > /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Counts_Gene_Activity/Adjusted/Gene_counts_"${Clusters[i]}"_Rep3.txt
#done

Controls=(Control_for_Cluster1 Control_for_Cluster2 Control_for_Cluster3 Control_for_Cluster4)
#for i in ${!Controls[@]}; do cut -f1,12 ${Controls[i]}_Rep1.metadata.txt > ${Controls[i]}_Rep1.barcodes.txt
#cut -f1,12 ${Controls[i]}_Rep2.metadata.txt > ${Controls[i]}_Rep2.barcodes.txt
#cut -f1,12 ${Controls[i]}_Rep3.metadata.txt > ${Controls[i]}_Rep3.barcodes.txt
#sinto filterbarcodes -b $BAM -c ${Controls[i]}_Rep1.barcodes.txt --outdir "$OUT" -p 12
#sinto filterbarcodes -b $BAM -c ${Controls[i]}_Rep2.barcodes.txt --outdir "$OUT" -p 12
#sinto filterbarcodes -b $BAM -c ${Controls[i]}_Rep3.barcodes.txt --outdir "$OUT" -p 12

#done

Controls=(Control_for_Cluster_1 Control_for_Cluster_2 Control_for_Cluster_3 Control_for_Cluster_4)

for i in ${!Controls[@]}; do
bedtools coverage -a "$GENES" -b "${Controls[i]}"rep1.bam > /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Counts_Gene_Activity/Adjusted/Gene_counts_"${Controls[i]}"_Rep1.txt
bedtools coverage -a "$GENES" -b "${Controls[i]}"rep2.bam > /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Counts_Gene_Activity/Adjusted/Gene_counts_"${Controls[i]}"_Rep2.txt
bedtools coverage -a "$GENES" -b "${Controls[i]}"rep3.bam > /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Counts_Gene_Activity/Adjusted/Gene_counts_"${Controls[i]}"_Rep3.txt
done
