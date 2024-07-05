#!/bin/bash -l

#SBATCH --ntasks=12
#SBATCH --mem-per-cpu=55G
#SBATCH --time=1-00:15:00
#SBATCH --job-name="009_Generate_control_regions.sh"
#SBATCH --output=std/009_Generate_control_regions.sh.out
#SBATCH -p seymourlab


module load bedtools
module load samtools

ACRs_INFO=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/ACRs_class.distance.txt
GENES=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/WN_genes.simple.bed

GROUP_ACRs=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level

CONTROL=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Identify_mappable_regions/Control_regions

GENIC_ACRs=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/Genic_ACRs.bed
PROXIMAL_ACRs=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/Proximal_ACRs.bed
DISTAL_ACRs=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/Distal_ACRs.bed

NONCODING_GENIC_ACRs=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/Noncoding_genic_ACRs.bed

CODING_GENIC_ACRs=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/Coding_genic_ACRs.bed

EXONS=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/exons.bed

chrom_info=/bigdata/seymourlab/idiaz026/REFERENCES/PWN_HAPA/WN_HAPA_modContig_genomeFile.txt


TMPDIR=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/Temp


ALL_WGS=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-07-04_WGS_variant_callingPWN/Alignments/Merged_bams/Correct_RG/Mark_dup/Collate_WGS.bam

# samtools sort $ALL_WGS -o $CONTROL/Collate_WGS.sort.bam -O bam -m 2G --threads 24
# samtools index $CONTROL/Collate_WGS.sort.bam --threads 24

# samtools view -b -q 30 $CONTROL/Collate_WGS.sort.bam --threads 24 > $CONTROL/Collate_WGS.sort.UM.bam
# #include regions with coverage > 20 
# bedtools genomecov -ibam $CONTROL/Collate_WGS.sort.UM.bam -split -bg | awk '$4 > 20' | cut -f1,2,3 > $CONTROL/WGS_coverage_20cov.bed

bedtools merge -i $CONTROL/WGS_coverage_20cov.bed | sed 's/_/-/' > "$CONTROL"/WGS_mappable_regions.bed


MAPPABLE_REGIONS="$CONTROL"/WGS_mappable_regions.bed

awk 'BEGIN { OFS= "\t"}; {print$1,$2-2000,$3+2000}' $GENES | awk 'BEGIN {OFS = "\t"}; { if ($2 >=1) {print $1,$2,$3 } }' > $CONTROL/genes.plus.2kb.flank.bed
awk 'BEGIN { OFS= "\t"}; {print$1,$2-2000,$2}' $GENES | awk 'BEGIN {OFS = "\t"}; { if ($2 >=1) {print $1,$2,$3 } }' > $CONTROL/upstream_2kb.bed
awk 'BEGIN { OFS= "\t"}; {print$1,$3,$3+2000}' $GENES > $CONTROL/downstream_2kb.bed


awk 'BEGIN { OFS= "\t"}; {print$1,$2,$3}' $MAPPABLE_REGIONS > $CONTROL/mappable_regions.bed


cat $CONTROL/genes.plus.2kb.flank.bed $DISTAL_ACR > $CONTROL/genes.plus2kb.flank.DISTAL_ACR.toexclude.bed
cat $CONTROL/upstream_2kb.bed $CONTROL/downstream_2kb.bed > $CONTROL/up2kb_down2kb.flank.bed

rm $CONTROL/upstream_2kb.bed
rm $CONTROL/downstream_2kb.bed

sort -k1,1 -k2,2n $CONTROL/genes.plus2kb.flank.DISTAL_ACR.toexclude.bed > $CONTROL/genes.plus2kb.flank.DISTAL_ACR.toexclude.sort.bed
sort -k1,1 -k2,2n $CONTROL/up2kb_down2kb.flank.bed > $CONTROL/up2kb_down2kb.flank.sort.bed


bedtools intersect -f 1.0 -e -a $CONTROL/up2kb_down2kb.flank.sort.bed -b $CONTROL/mappable_regions.bed | sort -k1,1 -k2,2n > $CONTROL/mappable_proximal_regions.bed

rm $CONTROL/genes.plus2kb.flank.DISTAL_ACR.toexclude.bed
rm $CONTROL/up2kb_down2kb.flank.bed

#build set of genic ACRs and exonic regions to exclude for noncoding genic acrs

cat $GENIC_ACRs $EXONS | cut -f1,2,3 | sort -k1,1 -k2,2n > "$CONTROL"/Exonic_regions_for_exclusion.bed



echo "create random set of control regions to compare to DISTAL_ACR"
bedtools shuffle -seed 123 -i $DISTAL_ACRs -g $chrom_info -excl $CONTROL/genes.plus2kb.flank.DISTAL_ACR.toexclude.sort.bed -incl $CONTROL/mappable_regions.bed > $CONTROL/distal_control_regions.bed

echo "create random set of control regions to compare to GENIC_ACR"
bedtools shuffle -seed 123 -i $GENIC_ACRs -g $chrom_info -excl $GENIC_ACRs -incl $GENES > $CONTROL/genic_control_regions.bed

bedtools shuffle -seed 123 -i $NONCODING_GENIC_ACRs -g $chrom_info -excl "$CONTROL"/Exonic_regions_for_exclusion.bed -incl $GENES > $CONTROL/noncoding_genic_control_regions.bed


bedtools shuffle -seed 123 -i $CODING_GENIC_ACRs -g $chrom_info -excl $CODING_GENIC_ACRs -incl $EXONS > $CONTROL/coding_genic_control_regions.bed


echo "create random set of control regions to compare to PROXIMAL_ACR"
bedtools shuffle -seed 123 -i $PROXIMAL_ACRs -g $chrom_info -excl $PROXIMAL_ACRs -incl $CONTROL/mappable_proximal_regions.bed > $CONTROL/proximal_control_regions.bed


#rm genes.plus2kb.flank.DISTAL_ACR.toexclude.sort.bed
#rm up2kb_down2kb.flank.sort.bed
#rm genes.plus.2kb.flank.bed
##########

cd $CONTROL
# echo "Further subset control regions"
class=(distal coding_genic noncoding_genic proximal)
FREQ=(Unique Common High_frequency)

for i in ${!class[@]}; do
      cat ${class[i]}_control_regions.bed | grep "${FREQ[0]}" > $CONTROL/${class[i]}_${FREQ[0]}_control.bed
      cat ${class[i]}_control_regions.bed | grep "${FREQ[1]}" > $CONTROL/${class[i]}_${FREQ[1]}_control.bed
      cat ${class[i]}_control_regions.bed | grep "${FREQ[2]}" > $CONTROL/${class[i]}_${FREQ[2]}_control.bed
done

cd $GROUP_ACRs

CLASS=(Distal Coding_genic Noncoding_genic Proximal)
for i in ${!CLASS[@]}; do
      cat ${CLASS[i]}_ACRs.bed | grep "${FREQ[0]}" > $GROUP_ACRs/${CLASS[i]}_${FREQ[0]}_ACRs.bed
      cat ${CLASS[i]}_ACRs.bed | grep "${FREQ[1]}" > $GROUP_ACRs/${CLASS[i]}_${FREQ[1]}_ACRs.bed
      cat ${CLASS[i]}_ACRs.bed | grep "${FREQ[2]}" > $GROUP_ACRs/${CLASS[i]}_${FREQ[2]}_ACRs.bed
done
