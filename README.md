# scATAC-Seq Analysis Pipeline

This repository provides a Snakemake pipeline for analyzing single-cell ATAC-seq (scATAC-seq) data using various tools, including `cellranger-atac`, `Socrates`, `Genrich`, and `bedtools`. The pipeline processes raw scATAC fastq files, performs barcode cleaning, demultiplexes cells, and identifies accessible chromatin regions (ACRs).

## Table of Contents
- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
- [Pipeline Steps](#pipeline-steps)
- [Outputs](#outputs)

## Installation
### Prerequisites
- Install [Snakemake](https://snakemake.readthedocs.io/)
- Required software tools: `cellranger-atac`, `sinto`, `samtools`, `R`, `popscle`, `bedtools`, and `ncbi-blast`
- Ensure access to SLURM for resource allocation on a high-performance computing (HPC) cluster.

### Set up Conda environments
For dependencies, create Conda environments as listed in `config/SocratesEnv.yaml` and `config/DemuxletEnv.yaml` 
conda env create -f config/SocratesEnv.yaml
conda env create -f config/DemuxletEnv.yaml

## Configuration
### Edit  paths in config/config.yaml
### for cellranger
1. CELLRANGER_PATH: path to cellranger-atac software
2. CELLRANGER_ref: path to reference created with cellranger
3. scATACraw: path to fastq of single-cell ATAC-seq reads
4. NAME: Sample name of reference for cellranger-atac
5. Nuclear: expression for scaffold prefix for nuclear genome for sinto usage
6. Plastid: expression for scaffold prefix for plastid genome
7. GFF: gene annotation in gff format
8. CHRFILE: lengths of scaffolds

### for demultiplexing
9. MACSpath: path to macs2 software
10. VCFdemux: path to VCF for demuxlet
11. SAMPLENAMES: file with atac-seq sample names
12. PICARD: path to picard software

### for ACR calling
13. GENRICH: path to genrich software
14. SAMPLEMETA: "Metadata/SampleNames.txt"
15. WGS: "/path/to/WGS_alignments"
16. REFERENCE: "/path/to/reference_genome.fa"
17. PLASTID: "/path/to/plastid_blast_db"


## Pipeline Steps

1. Cell Ranger ATAC Processing
Description: Runs cellranger-atac to process raw scATAC fastq files and generate BAM and fragment files.
Input: Raw fastq files, reference genome.
Output: {NAME}/outs/possorted_bam.bam, {NAME}/outs/fragments.tsv
2. Cell Barcode Cleaning
Description: Formats and filters cell barcodes using the sinto tool.
Input: BAM and fragment files from Cell Ranger ATAC output.
Output: Filtered fragments and contig lists.
3. Cell Identification (Socrates)
Description: Uses Socrates to filter cell barcodes based on genomic regions.
Output: Filtered barcodes and sparse matrices for further analysis.
4. Demultiplexing Cells
Description: Uses demuxlet with SNP data to demultiplex cells.
Input: BAM file, VCF file of SNPs, sample names.
Output: Best match demultiplex results, clean barcode list.
5. Filtering Barcodes
Description: Filters barcodes using picard to remove duplicates and produce a final BAM.
Output: Filtered and indexed BAM files.
6. Downsample Nuclei
Description: Randomly downsamples barcodes to a specified count (e.g., 552 per sample).
Input: List of clean barcodes.
Output: Downsampled barcode files.
7. ACR Calling
Description: Uses Genrich to call Accessible Chromatin Regions (ACRs) based on downsampled BAM files.
Output: Peak files and ACR log files.
8. ACR Processing
Description: Filters out regions aligning to plastid genomes and calculates jaccard similarity.
Output: Processed and concatenated ACR files for frequency analysis.
9. ACR Analysis and Visualization
Description: Analyzes and visualizes ACR positions relative to genes.
Output: Merged and classified ACR files, visualization plots.
10. ACR Classification
Description: Classifies ACRs based on proximity to genes and exons.
Output: ACRs classified by gene location, overlap with exons, and nearest gene.

## Outputs

The pipeline produces a variety of output files, including:

Filtered BAM files: Located in DEMUX/BAMscATAC/
Downsampled ACR files: Stored in DownsampledBams/Peaks/
Merged ACRs: Final classified ACR files for downstream analysis
Visualization files: Plots summarizing ACRs
ACRs classified by genomic context ACRs_classified/
