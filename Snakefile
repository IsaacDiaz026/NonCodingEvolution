configfile: "config/config.yaml"

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
        slurm_extra="--output=Logs/Socrates_iscell.out"
    conda:
        "config/SocratesEnv.yaml"
    params:
    name=config["NAME"] 
    demux=config["DEMUX"]
    input:
    frags="All_ATAC_plusplastidfragments.bed.gz"
    gff=config["GFF"] #annotation for Socrates
    chr=config["CHRFILE"] #contig lengths for Socrates
    macs=config["MACSpath"] #path to MACS2 software for ACR calling
    plastidcontig="plastid_contigs.txt"
    reads="{params.name}/outs/possorted_bam.bam"
    output:
    "Tn5_frag_conv2bed.bed",
    "All_ATAC.sparse",
    "raw_cell_data_isCELL.txt",
    "barcodes_to_keep.txt",
    "scATAC_filtered.bam"
    shell:
    """
    Rscript 001_socrates_isCell.R {input.frags} {input.gff} {input.chr} {input.macs} {input.plastidcontig} 
    gzip {output[2]} 
    cut -f1,18 {output[2]} | tail -n +2 > {output[3]}
    sinto filterbarcodes -b {input.reads} -c {output[3]} -p 16 --outdir {params.demux}
    mv "{params.demux}/A.bam" "{params.demux}/scATAC_filtered.bam" 
    """


