scATAC-Seq Analysis Pipeline

This repository provides a Snakemake pipeline for analyzing single-cell ATAC-seq (scATAC-seq) data using various tools, including cellranger-atac, Socrates, Genrich, and bedtools. The pipeline processes raw scATAC fastq files, performs barcode cleaning, demultiplexes cells, and identifies accessible chromatin regions (ACRs).

Table of Contents

Installation
Configuration
Usage
Pipeline Steps
Outputs
Installation

Prerequisites
Install Snakemake
Required software tools: cellranger-atac, sinto, samtools, R, popscle, bedtools, and ncbi-blast
Ensure access to SLURM for resource allocation on a high-performance computing (HPC) cluster.
Set up Conda environments
For dependencies, create Conda environments as listed in config/SocratesEnv.yaml and config/DemuxletEnv.yaml. Example:


conda env create -f config/SocratesEnv.yaml
conda env create -f config/DemuxletEnv.yaml
Configuration

The pipeline uses a configuration file located at config/config.yaml. Modify the paths to suit your data and environment.

Example config.yaml structure:

# cell-ranger atac usage
CELLRANGER_PATH: "/path/to/cellranger-atac"
CELLRANGER_ref: "/path/to/reference"
scATACraw: "/path/to/raw_scATAC_fastq"
NAME: "SampleName"

# Sinto barcode cleaning
Nuclear: "?i)^Scaffold_"
Plastid: "?i)^NC-"

# Socrates and ACR calling
GFF: "/path/to/genes.gff"
CHRFILE: "/path/to/contig_lengths.txt"
VCFdemux: "/path/to/vcf_file.vcf.gz"
PICARD: "/path/to/picard.jar"

# ACR calling and analysis
GENRICH: "/path/to/Genrich"
SAMPLEMETA: "MetadataSampleNames.txt"
WGS: "/path/to/WGS_alignments"
REFERENCE: "/path/to/reference_genome.fa"
PLASTID: "/path/to/plastid_db"


Create the following conda environments using yaml files:
snakemake --configfile config/config.yaml --cores <number_of_cores>



Pipeline Steps

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
Outputs

The pipeline produces a variety of output files, including:

Filtered BAM files: Located in DEMUX/BAMscATAC/
Downsampled ACR files: Stored in DownsampledBams/Peaks/
Merged ACRs: Final classified ACR files for downstream analysis
Visualization files: Plots summarizing ACRs
This documentation should provide an accessible overview of each step, with paths configurable as needed for different systems. Let me know if you'd like additional customization!