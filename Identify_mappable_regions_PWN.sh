#!/bin/bash -l

#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=300G
#SBATCH --job-name="Identification_of_Mappable_regions_Cclem.sh"
#SBATCH --output=std/Identification_of_Mappable_regions_Cclem.sh.out
#SBATCH -p highmem

#This script is from Lu et al. 2019 "Prevalence evolution and chromatin signatures of plant regulatory elements"

module load bedtools
module load bowtie
module load bwa
module load samtools

chrom_info=/bigdata/seymourlab/idiaz026/REFERENCES/Cclementina/Cclementina_genomeFile.txt
fasta=/bigdata/seymourlab/idiaz026/REFERENCES/Cclementina/Cclementina_182_v1.fa
len=150
INDEX=/bigdata/seymourlab/idiaz026/REFERENCES/Cclementina/Cclementina_182_v1
thread=1
one=100000000      ##### max limit of one reads set
win=1000000        ##### max limit to merge

GENES=/bigdata/seymourlab/idiaz026/REFERENCES/Cclementina/Cclem_longest_mRNA.sort.gff3



cd /bigdata/seymourlab/idiaz026/Results/Floral_scATAC/Identify_mappable_regions

echo -n > SE$len.bed
cat $chrom_info | while read LINE
do
      chr=$(echo "$LINE"| cut -f1)
      end=$(echo "$LINE"| cut -f2)
      seq $len $end > SE$len.tmp000
      echo -n >SE$len.$chr.sim.bed
      awk '{print "'$chr'",$1-'$len',$1}' OFS="\t" SE$len.tmp000 >> SE$len.$chr.sim.bed
      rm SE$len.tmp000
      m=$(($end/$one))
      split -l $one -d -a ${#m} SE$len.$chr.sim.bed SE$len.$chr.sim.bed_split_
      for n in `seq -w 0 $m`
      do
            bedtools getfasta -fi $fasta -bed SE$len.$chr.sim.bed_split_$n -fo SE$len.$chr.sim.bed_split_$n.fa
            bowtie $INDEX SE$len.$chr.sim.bed_split_$n.fa -S SE$len.$chr.split_$n.sam -f -t -p $thread -v 2 --best --strata -m 1
            samtools sort -O 'bam' -o SE$len.$chr.split_$n.bam -T tmp$chr.split_$n SE$len.$chr.split_$n.sam
            rm SE$len.$chr.split_$n.sam SE$len.$chr.sim.bed_split_$n.fa SE$len.$chr.sim.bed_split_$n
            bedtools bamtobed -i SE$len.$chr.split_$n.bam >SE$len.$chr.split_$n.bed
            rm SE$len.$chr.split_$n.bam
      done
      cat SE$len.$chr.split_*.bed > SE$len.$chr.bed
      rm SE$len.$chr.split_*.bed SE$len.$chr.sim.bed
      cat SE$len.$chr.bed >> SE$len.bed

      l=$(cat SE$len.$chr.bed|wc -l)
      x=$(($l/$win))
      split -l $win -d -a ${#x} SE$len.$chr.bed SE$len.$chr.bed_split_

      for y in `seq -w 0 $x`
      do
      bedtools merge -i SE$len.$chr.bed_split_$y > SE$len.$chr.bed_$y.x
      rm SE$len.$chr.bed_split_$y
      done

      cat SE$len.$chr.bed_*.x |bedtools merge -i - >  SE$len.$chr.mappable.bed
      rm SE$len.$chr.bed_*.x
      done

      bedtools sort -i SE$len.bed | bedtools genomecov -i - -split -bg -g $chrom_info > SE$len.bg
