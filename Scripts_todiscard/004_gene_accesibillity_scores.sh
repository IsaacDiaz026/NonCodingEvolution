#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --mem=150G
#SBATCH --time=10-00:15:00
#SBATCH --job-name="004_gene_accesibillity_scores.sh"
#SBATCH --output=std/004_gene_accesibillity_scores.sh.out
#SBATCH -p highmem

conda activate Socrates

Rscript gene_accessibillity_scores.R
