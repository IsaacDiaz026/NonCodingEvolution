#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=1-00:15:00
#SBATCH --job-name="002b_peak_calling.sh"
#SBATCH --output=std/002b_peak_calling.sh.sh.out
#SBATCH -p seymourlab 

conda activate Socrates_fix
module load samtools

module load bedtools
module load ncbi-blast/2.2.30+

BAMS=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/SintoGenotypeBams/Collate_BAM
TMPDIR=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Temp

path_MtPT_Db=/bigdata/seymourlab/idiaz026/Fairchild/Assembly/Release1.0/References/Plastid_mit_db

GENOME_FILE=/bigdata/seymourlab/idiaz026/REFERENCES/PWN_HAPA/WN_HAPA_modContig_genomeFile.txt
REFERENCE=/bigdata/seymourlab/idiaz026/REFERENCES/PWN_HAPA/WN_HapA_modContig.fa

cd /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/SintoGenotypeBams/

#samtools index "$BAMS"/All_ATAC.clean.bam

#sinto fragments -b "$BAMS"/All_ATAC.clean.bam -f "$BAMS"/All_ATAC_fragments.bed --use_chrom "(?i)^Scaffold-" --collapse_within

#gzip "$BAMS"/All_ATAC_fragments.bed

#Rscript /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/002b_peakcalling.R

#gzip /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/All_ATAC_peaks.CLEAN.sparse

PEAKS=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_by_geno/Downsampled_genrich

cd "$PEAKS"

# for file in *.txt.bed ; do sort -k1,1 -k2,2n "$file" > "$file".sort.bed ; bedtools getfasta -fi "$REFERENCE" -bed "$file".sort.bed  > "$file".fa ; blastn -db $path_MtPT_Db -query "$file".fa -outfmt 7 -out "$file"_regions.mtpt ; sed '/#/d' "$file"_regions.mtpt | cut -f1 | sed "s/\:/\t/g" | sed "s/-/\t/2" | bedtools sort -i - | uniq > "$file"_open_regions.black ; bedtools intersect -a "$file".sort.bed -b "$file"_open_regions.black  -v | cut -f1-4 > "$file".sort.CLEAN.bed ; done


# parallel "bedtools jaccard -a {1} -b {2} -g /bigdata/seymourlab/idiaz026/REFERENCES/PWN_HAPA/WN_HAPA_modContig_genomeFile.txt | awk 'NR>1' | cut -f 3 > {1}.{2}.jaccard" ::: `ls *.sort.CLEAN.bed` ::: `ls *.sort.CLEAN.bed`

# find . \
# | grep jaccard \
# | xargs grep "" \
# | sed -e s"/\.\///" \
# | perl -pi -e "s/.sort.CLEAN.bed./.sort.CLEAN.bed\t/" \
# | perl -pi -e "s/.jaccard:/\t/" \
# > pairwise.ATAC.txt

# cat pairwise.ATAC.txt \
# | sed -e 's/.05.peaks.txt.bed//g' \
# > pairwise.ATAC.shortnames.txt


#Rscript /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Jaccard_similarity.R

OVERLAP=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACR_OVERLAP

# bedtools multiinter -i "$PEAKS"/*.sort.CLEAN.bed \
#      | awk '{print $4"\t"$3-$2}' \
#      | sort -k1,1n \
#      | bedtools groupby -g 1 -c 2 -o sum \
# > "$OVERLAP"/atac.occupancy.dist.txt

# cd "$OVERLAP"

# cat "$PEAKS"/*.sort.CLEAN.bed > "$OVERLAP"/Cat_all_acrs.bed

# cd $PEAKS 

# cat SCFS_citron*.sort.CLEAN.bed Microcitrus_australiasica*.sort.CLEAN.bed C_macrophylla*.sort.CLEAN.bed Algerian_clementine*.sort.CLEAN.bed Kao_panne_pummelo*.sort.CLEAN.bed > "$OVERLAP"/Cat_species_acrs.bed


# cd "$OVERLAP"

# sort -k1,1 -k2,2n Cat_species_acrs.bed > Cat_species_acrs.sort.bed

# sort -k1,1 -k2,2n Cat_all_acrs.bed > Cat_all_acrs.sort.bed

# bedtools merge -i Cat_all_acrs.sort.bed > Merged_ACRs.bed

# bedtools merge -i Cat_species_acrs.sort.bed > Merged_species_ACRs.bed

# bedtools intersect -wa -C -names Algerian_clementine Bouquettier_de_Nice Citrus_junos C_macrophylla Frost_lisbon Frost_owen_satsuma Kao_panne_pummelo Koster_mandarin Microcitrus_australiasica Orlando Oroblanco_grpft SCFS_citron Seedless_kishu S.emperor_ponkan Swingle Tango USDA_88-2 Valentine Washington_navel -a Merged_ACRs.bed -b "$PEAKS"/*.sort.CLEAN.bed > Merged_ACRs.counts.bed

# bedtools intersect -wa -C -names SCFS_citron Microcitrus_australiasica C_macrophylla Algerian_clementine Kao_panne_pummelo -a Merged_species_ACRs.bed -b "$PEAKS"/SCFS_citron*.sort.bed "$PEAKS"/Microcitrus_australiasica*.sort.CLEAN.bed "$PEAKS"/C_macrophylla*.sort.CLEAN.bed "$PEAKS"/Algerian_clementine*.sort.CLEAN.bed "$PEAKS"/Kao_panne_pummelo*.sort.CLEAN.bed > Merged_ACRs.counts.species.bed


# Rscript /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/006a_plot_ACR_position_overlap.R

# cat /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/WN_HAPA_genes_modCTG.bed | cut -f1,2,3 > /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/WN_genes.simple.bed
GENES=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/WN_genes.simple.bed

cd /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/

sort -k1,1 -k2,2n /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/ACRs_classified.bed > /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/ACRs_classified.sort.bed


bedtools closest -d -a /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/ACRs_classified.sort.bed -b "$GENES" > /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/ACRs_class.distance.txt

#species level
sort -k1,1 -k2,2n /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/ACRs_classified.bed > /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/ACRs_classified.sort.bed

bedtools closest -d -a /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/ACRs_classified.sort.bed -b "$GENES" > /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/ACRs_class.distance.txt

awk 'BEGIN {OFS="\t"}; {if ($8 == 0) {print $1,$2,$3,$4} }' /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/ACRs_class.distance.txt > /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/genic_to_specify.bed

awk 'BEGIN {OFS="\t"}; {if ($8 == 0) {print $1,$2,$3,$4} }' /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/ACRs_class.distance.txt > /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/genic_to_specify.bed

cat /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/WN_HAPA_genes_modCTG.gff | grep "CDS" | sort -k1,1 -k4,4n | bedtools merge -i - > /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/exons.bed 

bedtools intersect -wao -a /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/genic_to_specify.bed -b /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/exons.bed > /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/genic_acrsIntWexons.txt

bedtools intersect -wao -a /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/genic_to_specify.bed -b /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/exons.bed > /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/genic_acrsIntWexons.txt


#Rscript /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/007_Describe_ACRs.R

echo 'get closest gene (mRNA only) to every ACR using including strand info'

GFF=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/WN_HAPA_genes_modCTG.gff

grep "mRNA" $GFF > /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/PWN_modCTG_mRNAonly.gff

MRNA=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/PWN_modCTG_mRNAonly.gff


ALL_GENO_ACRS=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/ACRs_classified.sort.bed
SPECIES_ACRS=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/ACRs_classified.sort.bed


bedtools closest -D b -a $ALL_GENO_ACRS -b "$MRNA" > /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/ACRs_closest.mRNA.txt

bedtools closest -D b -a $SPECIES_ACRS -b "$MRNA" > /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/ACRs_closest.mRNA.txt




