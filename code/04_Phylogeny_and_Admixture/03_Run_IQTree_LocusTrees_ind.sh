#!/bin/bash

# 02_Run_Prep_For_IQTree.sh
# Dec 2023
# JRG
# Submits to queue for Phylogeny Reconstruction

root="$(git rev-parse --show-toplevel)"
source "${root}/pkings_gwas.env"

mkdir -p "${root}/output_data/slurm_logs/04_Phylogeny_And_Admixture/"

# Default array jobs range
default_range="0-500"

# If the first argument is not empty, use it as the array jobs range
# Otherwise, use the default range
array_jobs=${1:-$default_range}

echo "sbatch --job-name \"PHYLO_IQTREE_LOCI\" --output ${root}/output_data/slurm_logs/04_Phylogeny_And_Admixture/phylo_run_iqtree_loci_extracted_%a.log -a ${array_jobs} --export=root=${root} ${root}/code/04_Phylogeny_and_Admixture/run_iqtree_gene_loci_ind.sb"
sbatch --job-name "PHYLO_IQTREE_LOCI" --output "${root}/output_data/slurm_logs/04_Phylogeny_And_Admixture/phylo_run_iqtree_loci_extracted_%a.log" -a "${array_jobs}" --export=root=${root} "${root}/code/04_Phylogeny_and_Admixture/run_iqtree_gene_loci_ind.sb"
