#!/bin/bash

# 02_Run_Prep_For_IQTree.sh
# Dec 2023
# JRG
# Submits to queue for Phylogeny Reconstruction

root="$(git rev-parse --show-toplevel)"
source "${root}/pkings_gwas.env"

mkdir -p "${root}/output_data/slurm_logs/04_Phylogenomics/"

echo "sbatch --job-name \"PHYLO_IQTREE_LOCI\" --output ${root}/output_data/slurm_logs/04_Phylogenomics/phylo_run_iqtree_loci_together.log --export=root=${root} ${root}/code/04_Phylogenomics/run_iqtree_gene_loci_together.sb"
sbatch --job-name "PHYLO_IQTREE_LOCI" --output "${root}/output_data/slurm_logs/04_Phylogenomics/phylo_run_iqtree_loci_together.log" --export=root=${root} "${root}/code/04_Phylogenomics/run_iqtree_gene_loci_together.sb"
