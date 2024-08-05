configfile: "config/config.yaml"

import pandas as pd

samples = pd.read_table("Metadata/SampleNames.txt").set_index("genotype", drop=False)
all_ids=list(samples["genotype"])


rule cellranger:
    """
    Run cellranger atac to process raw scATAC fastq files
    must run mkref to produce cellranger compatible reference
    """
    resources:
        slurm_partition="highmem",
        tasks=32,
        mem_mb=700000,
        nodes=1,
        runtime=200000
        slurm_extra="--output=Logs/cellranger.out"
    input:
        ref=config["CELLRANGER_ref"], # path to reference directory from mkref
        reads=config["scATACraw"]
    params:
        cellrangerpath=config["CELLRANGER_PATH"],
        name=config["NAME"]
    log:
        "logs/cellrangerCount.log"
    output:
        "{params.name}/outs/possorted_bam.bam",
        "{params.name}/outs/fragments.tsv"
    shell: 
        "{params.cellrangerpath} count --id={params.name} --reference={input.ref} --fastqs={input.reads}"


rule CleanCells:
    """"
    Format scATAC data for R package Socrates, then filter the scATAC cells
    """"
    resources:
        slurm_partition="highmem",
        tasks=16,
        mem_mb=900000,
        nodes=1,
        runtime=200000,
        slurm_extra="--output=Logs/cleancells.out"
    conda:
        "config/SocratesEnv.yaml"
    params:
        Nuclear=config["Nuclear"],
        Plastid=config["Plastid"],
        name=config["NAME"]  
    input:
        reads="{params.name}/outs/possorted_bam.bam",
        frags="{params.name}/outs/fragments.tsv"
    output:
        "All_ATAC_plastidfragments.bed",
        "plastid_contigs.txt",
        "All_ATAC_plusplastidfragments.bed.gz"
    shell:
        """
        sinto fragments -b {input.reads} -f {output[0]} -p 16 --use_chrom "{params.plastid}" 
        cut -f1 {output[0]} | uniq > {output[1]}
        cat {input.frags} {output[0]} | gzip > {output[2]}
        """


rule IdentifyCells:
    """
    Use socrates to filter cell barcodes
    """
    resources:
        slurm_partition="highmem",
        tasks=16,
        mem_mb=900000,
        nodes=1,
        runtime=200000,
        slurm_extra="--output=Logs/IdentifyCells.out"
    conda:
        "config/SocratesEnv.yaml"
    params:
        name=config["NAME"] .
        demux=config["DEMUX"]
    input:
        frags="All_ATAC_plusplastidfragments.bed.gz",
        gff=config["GFF"], #annotation for Socrates
        chr=config["CHRFILE"], #contig lengths for Socrates
        macs=config["MACSpath"], #path to MACS2 software for ACR calling
        plastidcontig="plastid_contigs.txt",
        reads="{params.name}/outs/possorted_bam.bam"
    output:
        "Tn5_frag_conv2bed.bed",
        "All_ATAC.sparse",
        "raw_cell_data_isCELL.txt",
        "DEMUX/barcodes_to_keep.txt",
        "DEMUX/scATAC_filtered.bam"
    shell:
        """
        TMPDIR=Temp
        mkdir -p DEMUX
        Rscript 001_socrates_isCell.R {input.frags} {input.gff} {input.chr} {input.macs} {input.plastidcontig} 
        gzip {output[2]} 
        cut -f1,18 {output[2]} | tail -n +2 > {output[3]}
        sinto filterbarcodes -b {input.reads} -c {output[3]} -p 16 --outdir DEMUX
        mv DEMUX/A.bam DEMUX/scATAC_filtered.bam
        """


rule Demultiplex:
    """
    Use demuxlet to demultiplex scATAC using SNPs
    Intersect with clean barcodes detected with Socrates
    """
    resources:
        slurm_partition="highmem",
        tasks=16,
        mem_mb=600000,
        nodes=1,
        runtime=200000,
        slurm_extra="--output=Logs/Demultiplex.out"
    conda:
        "config/DemuxletEnv.yaml"
    params:
        name=config["NAME"] 
    input:
        reads="{params.name}/outs/possorted_bam.bam",
        variants=config["VCF"],
        barcodes="DEMUX/barcodes_to_keep.txt",
        samples=config["SAMPLENAMES"]
    output:
        "{params.name}".best,
        "DEMUX/barcodesClean_minusInter.txt"

    shell:
        """
        TMPDIR=Temp
        popscle demuxlet --sam {input.reads} --vcf {input.variants} --field GT --out {output[0]}
        Rscript 002_Demuxlet_eval.R {output[0]} {input.barcodes} {input.samples}
        """


rule FilterBarcodes:
    """
    Run demultiplexing using clean barcodes
    """
    resources:
        slurm_partition="highmem",
        tasks=16,
        mem_mb=600000,
        nodes=1,
        runtime=200000,
        slurm_extra="--output=Logs/FiltBarcodes.out"
    conda:
        "config/SocratesEnv.yaml"
    params:
        name=config["NAME"] 
        picard=config["PICARD"]
        Nuclear=config["Nuclear"]
    input:
        reads="{params.name}/outs/possorted_bam.bam",
        variants=config["VCF"],
        barcodes="DEMUX/barcodesClean_minusInter.txt",
        samples=config["SAMPLENAMES"]
    output:
        "DEMUX/BAMscATAC/Tango.bam",
        "DEMUX/BAMscATAC/All_ATAC.bam",
        "DEMUX/BAMscATAC/All_ATAC.rmdup.bam",
        "DEMUX/BAMscATAC/All_ATAC_rmdup_metrics.txt",
        "DEMUX/BAMscATAC/All_ATAC.clean.bam",
        "DEMUX/BAMscATAC/All_ATAC_fragments.bed"
    shell:
        """
        TMPDIR=Temp
        mkdir -p DEMUX/BAMscATAC

        sinto filterbarcodes -b {input.reads} -p 16 -c {input.barcodes} --outdir DEMUX/BAMscATAC
        samtools index DEMUX/BAMscATAC/*bam

        samtools merge {output[1]} DEMUX/BAMscATAC/*bam
        samtools index {output[1]}

        java -jar {params.picard} MarkDuplicates INPUT= {output[1]} \
        O= {output[2]} M= {output[3]} \
        REMOVE_DUPLICATES=true ASSUME_SORTED=false \
        VALIDATION_STRINGENCY=LENIENT BARCODE_TAG=CB

        samtools index {output[2]}

        samtools view -f 3 -bhq 10 {output[2]} > {output[4]}

        samtools index {output[4]}

        sinto fragments -b {output[4]} -f {output[5]} -p 16 --use_chrom {params.Nuclear} --collapse_within

        gzip {output[5]}
        mkdir PeakCalling

        """



rule DownsampleNuclei:
    """
     "If number of barcodes is greater than 552, 
     downsample 552 barcodes randomly. 
     552 is median number of barcodes across 19 genotypes"
    """
    resources:
        slurm_partition="intel",
        tasks=1,
        mem_mb=10000,
        nodes=1,
        runtime=200000,
        slurm_extra="--output=Logs/downsampleNuc.out"
    conda:
        "config/SocratesEnv.yaml"
    threads: 20
    input:
        geno=config[SAMPLEMETA],
        cells="DEMUX/barcodesClean_minusInter.txt"
    expand:
    output:
        Barcodes/{sample}.barcodes.txt,
        Barcodes/{sample}.ds_barcodes.txt
    shell:
        """
        TMPDIR=Temp
        cut -f1 {input.geno} | tail -n + 2 > samples

        while read -r line; do 
            cat {input.cells} | grep "$line" > Barcodes/"$line".barcodes.txt
            num_barcodes=$(cat "$line".barcodes.txt | wc -l)
            if [["$num_barcodes" -gt 552]]; then
                 cat Barcodes/"$line".barcodes.txt | shuf -n 552 > Barcodes/"$line".ds_barcodes.txt
            else
                cat Barcodes/"$line".barcodes.txt > Barcodes/"$line".ds_barcodes.txt
            fi
        done < samples
        """

rule ACRcalling:
    """
    extract alignments from downsampled barcodes
    call ACRs with genrich
    """
    resources:
        slurm_partition="intel",
        tasks=20,
        mem_mb=80000,
        nodes=1,
        runtime=200000,
        slurm_extra="--output=Logs/Peakcall.out"
    conda:
        "config/SocratesEnv.yaml"
    params:
        name=config["NAME"], 
        picard=config["PICARD"],
        Nuclear=config["Nuclear"],
        genrich=config["GENRICH"],
    threads: 20
    input:
        reads="DEMUX/BAMscATAC/All_ATAC.clean.bam",
        barcodes="Barcodes/{sample}.ds_barcodes.txt",
        wgscontrol=config[WGS],
        samples=config["SAMPLENAMES"]
    output:
        "DownsampledBams/{sample}.bam",
        "$TMPDIR/{sample}.ATAC.namesort.bam",
        "$TMPDIR/{sample}.wgs.namesort.bam",
        "DownsampledBams/Peaks/{sample}.05.peaks.txt"
    shell:
        """
        TMPDIR=Temp
        mkdir -p PCR_DUP
        mkdir -p DownsampledBams
        mkdir -p DownsampledBams/Peaks
        sinto filterbarcodes -b {input.reads} -c {input.barcodes} \
        -p 8 --outdir DownsampledBams 
        
        samtools sort -n {ouput[0]} -o {output[1]}
        samtools sort -n "{input.wgscontrol}/{sample}.marked_duplicates.bam" -o {output[2]}

        {params.genrich} -t {output[1]} -R PCR_DUP/{sample}.pcrdups.txt -p 0.05 -r -y -j -o {output[3]} -f DownsampledBams/Peaks/{sample}.05.signal.log -c {output[2]}

        """

rule ACRprocessing:
    """
    remove ACRs that align to plastid genomes
    Classify ACRs and concatenate them for frequency analysis
    Need bedtools, ncbiblast, parallel
    """
    resources:
        slurm_partition="intel",
        tasks=1,
        mem_mb=50000,
        nodes=1,
        runtime=200000,
        slurm_extra="--output=Logs/ACRprocess.out"
    conda:
        "config/SocratesEnv.yaml"
    threads: 20
    input:
        samples=config["SAMPLENAMES"],
        reference=config["REFERENCE"],
        plastid=config["PLASTID"],
        genomefile=config["CHRFILE"]
    output:
        "pairwise.ATAC.shortnames.txt",
        "ACR_OVERLAP/atac.occupancy.dist.txt",
        "ACR_overlap/Cat_all_acrs.bed",
        "ACR_overlap/Cat_species_acrs.bed",
    shell:
        """
        module load bedtools
        module load ncbi-blast/2.2.30+
        module load parallel

        TMPDIR=Temp
        mkdir -p ACR_OVERLAP

        cd DownsampledBams/Peaks

        echo "Filtering ACRs that have chloroplast or mitochondria alignment"
        for file in *peaks.txt ; do sort -k1,1 -k2,2n "$file" | cut -f1,2,3,4 > "$file".sort.bed
        bedtools getfasta -fi {input.reference} -bed "$file".sort.bed > "$file".fa
        blastn -db {input.plastid} -query "$file".fa -outfmt 7 -out "$file"_regions.mtpt
        sed '/#/d' "$file"_regions.mtpt | cut -f1 | sed "s/\:/\t/g" | sed "s/-/\t/2" \
        | bedtools sort -i - | uniq > "$file"_open_regions.black
        bedtools intersect -a "$file".sort.bed -b "$file"_open_regions.black -v | cut -f1-4 > "$file".sort.CLEAN.bed ; done
        
        #jaccard similarity
        parallel "bedtools jaccard -a {1} -b {2} \
        -g {input.genomefile} | awk 'NR>1' | cut -f 3 > {1}.{2}.jaccard" \
        ::: `ls *.sort.CLEAN.bed` ::: `ls *.sort.CLEAN.bed`

        find . \
        grep jaccard \
        | xargs grep "" \
        | sed -e s"/\.\///" \
        | perl -pi -e "s/.sort.CLEAN.bed./.sort.CLEAN.bed\t/" \
        | perl -pi -e "s/.jaccard:/\t/" \
        > pairwise.ATAC.txt

        cat pairwise.ATAC.txt \
        | sed -e 's/.05.peaks.txt//g' \
        > {output[0]}
        
        
        bedtools multiinter -i *.sort.CLEAN.bed \
        | awk '{print $4"\t"$3-$2}' \
        | sort -k1,1n \
        | bedtools groupby -g 1 -c 2 -o sum \
        > {output[1]}

        echo 'combine all 19 genotypes'
        cat *.sort.CLEAN.bed > {output[2]}

        
        echo 'combine species representatives'
        cat SCFS_citron*.sort.CLEAN.bed Microcitrus_australiasica*.sort.CLEAN.bed C_macrophylla*.sort.CLEAN.bed Algerian_clementine*.sort.CLEAN.bed Kao_panne_pummelo*.sort.CLEAN.bed > {output[3]}

        """


rule ACRanalysis:
    """
    Launch ACR analyses and visualization
    Need bedtools
    """
    resources:
        slurm_partition="intel",
        tasks=1,
        mem_mb=50000,
        nodes=1,
        runtime=200000,
        slurm_extra="--output=Logs/ACRanalysis.out"
    conda:
        "config/SocratesEnv.yaml"
    params:
        Allsamples="Total"
        Specieslevel="Species_level"
    threads: 20
    input:
        ACRs="Cat_all_acrs.bed",
        speciesACRs="Cat_species_acrs.bed",
        jaccard="pairwise.ATAC.shortnames.txt",
        occupancy="atac.occupancy.dist.txt"
    output:
        "Merged_ACRs.bed",
        "Merged_species_ACRs.bed",
        "Merged_ACRs.counts.bed",
        "Merged_ACRs.counts.species.bed" 
    shell:
        """
        module load bedtools
  
        TMPDIR=Temp

        mkdir -p ACRs_classified

        cd DownsampledBams/Peaks/ACR_overlap

        sort -k1,1 -k2,2n Cat_all_acrs.bed > Cat_all_acrs.sort.bed
        sort -k1,1 -k2,2n Cat_species_acrs.bed > Cat_species_acrs.sort.bed
        bedtools merge -i Cat_all_acrs.sort.bed > {output[0]}
        bedtools merge -i Cat_species_acrs.sort.bed > {output[1]}

        #all accesions
        bedtools intersect -wa -C -names Algerian_clementine Bouquettier_de_Nice Citrus_junos \
        C_macrophylla Frost_lisbon Frost_owen_satsuma Kao_panne_pummelo Koster_mandarin \
        Microcitrus_australiasica Orlando Oroblanco_grpft SCFS_citron Seedless_kishu \
        S.emperor_ponkan Swingle Tango USDA_88-2 Valentine Washington_navel \
        -a {output[0]} -b ../*.sort.CLEAN.bed > {output[2]}
        
        #Species level
        bedtools intersect -wa -C -names SCFS_citron Microcitrus_australiasica C_macrophylla \
        Algerian_clementine Kao_panne_pummelo -a Merged_species_ACRs.bed \
        -b ../SCFS_citron*.sort.bed ../Microcitrus_australiasica*.sort.CLEAN.bed \
        ../C_macrophylla*.sort.CLEAN.bed ../Algerian_clementine*.sort.CLEAN.bed \
        ../Kao_panne_pummelo*.sort.CLEAN.bed > {output[3]}

        Rscript ../../003_plot_ACR_position_overlap.R {input.jaccard} \
        {input.occupancy} {output[2]} {input.ACRs} {params.Allsamples}

        Rscript ../../003_plot_ACR_position_overlap.R {input.jaccard} \
        {input.occupancy} {output[3]} {input.speciesACRs} {params.Species_level}

        """


rule ACRclassification:
    """
    Classify species level ACRs based on their location relative to genes
    Need bedtools
    """
    resources:
        slurm_partition="intel",
        tasks=1,
        mem_mb=50000,
        nodes=1,
        runtime=200000,
        slurm_extra="--output=Logs/ACRclassification.out"
    conda:
        "config/SocratesEnv.yaml"
    params:
        Allsamples="Total"
        Specieslevel="Species_level"
    threads: 20
    input:
        GENES=config["GFF"],
        consensusACRs="Species_level_ACRs_classified.bed",
        allACRs="Species_level_CAT_ALLACRs_classified.bed"
        jaccard="pairwise.ATAC.shortnames.txt",
        occupancy="atac.occupancy.dist.txt"
    output:
        "ACRs_classified/ACRs_classified.sort.bed",
        "ACRs_classified/ACRs_class.distance.txt",
        "ACRs_classified/genic_to_specify.bed",
        "ACRs_classified/genic_acrsIntWexons.txt",
        "ACRs_classified/ACRs_closest.mRNA.txt"
        
    shell:
        """
        TMPDIR=Temp
        cat {input.GENES} | cut -f1,2,3 > WN_genes.simple.bed
        #sort ACRs
        sort -k1,1 -k2,2n {input.ACRs} > {output[0]}

        #get ACRs distance to nearest gene
        bedtools closest -d -a ACRs_classified/ACRs_classified.sort.bed -b WN_genes.simple.bed > {output[1]}
        
        #get genic ACRs
        awk 'BEGIN {OFS="\t"}; {if ($8 == 0) {print $1,$2,$3,$4} }' {output[1]} > {output[2]}

        cat {input.GENES} | grep "CDS" | sort -k1,1 -k4,4n | bedtools merge -i - > exons.bed

        #get ACRs overlap with CDS
        bedtools intersect -wao -a exons.bed > {output[3]}


        echo 'get closest gene (mRNA only) to every ACR using including strand info'

        grep "mRNA" {input.GENES} > PWN_modCTG_mRNAonly.gff

        bedtools closest -D b -a {output[0]} -b PWN_modCTG_mRNAonly.gff  > {output[4]}

        Rscript 007_Describe_ACRs.R
        """



#Rscript /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/007_Describe_ACRs.R

echo 'get closest gene (mRNA only) to every ACR using including strand info'

GFF=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/WN_HAPA_genes_modCTG.gff

grep "mRNA" $GFF > /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/PWN_modCTG_mRNAonly.gff

MRNA=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/PWN_modCTG_mRNAonly.gff


ALL_GENO_ACRS=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/ACRs_classified.sort.bed
SPECIES_ACRS=/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/ACRs_classified.sort.bed


bedtools closest -D b -a $ALL_GENO_ACRS -b "$MRNA" > /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/ACRs_closest.mRNA.txt

bedtools closest -D b -a $SPECIES_ACRS -b "$MRNA" > /bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/ACRs_closest.mRNA.txt
