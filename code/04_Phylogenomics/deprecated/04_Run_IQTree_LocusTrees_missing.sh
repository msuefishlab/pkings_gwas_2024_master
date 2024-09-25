#!/bin/bash

# 02_Run_Prep_For_IQTree.sh
# Dec 2023
# JRG
# Submits to queue for Phylogeny Reconstruction

root="$(git rev-parse --show-toplevel)"
source "${root}/pkings_gwas.env"

mkdir -p "${root}/output_data/slurm_logs/04_Phylogeny_And_Admixture/missing_genes"

echo "sbatch --job-name \"PHYLO_IQTREE_LOCI\" --output ${root}/output_data/slurm_logs/04_Phylogeny_And_Admixture/missing_genes/phylo_run_iqtree_loci_%a.log -a 3 --export=root=${root} ${root}/code/04_Phylogeny_and_Admixture/run_iqtree_gene_loci_ind_missing.sb"
sbatch --job-name "PHYLO_IQTREE_LOCI" --output "${root}/output_data/slurm_logs/04_Phylogeny_And_Admixture/missing_genes/phylo_run_iqtree_loci_%a.log" -a 3 --export=root=${root} "${root}/code/04_Phylogeny_and_Admixture/run_iqtree_gene_loci_ind_missing.sb"
