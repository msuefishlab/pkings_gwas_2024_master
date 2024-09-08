#!/bin/bash

# 02_Run_Prep_For_IQTree.sh
# Dec 2023
# JRG
# Submits to queue for Phylogeny Reconstruction

root="$(git rev-parse --show-toplevel)"
source "${root}/pkings_gwas.env"

mkdir -p "${root}/output_data/slurm_logs/04_Phylogeny_And_Admixture/"

alignments=("${root}/output_data/04_Phylogeny_And_Admixture/1000_rand_filtered_gene_alignments/"*.fasta)
nalignments=${#alignments[@]}

# Default array jobs range
default_range="0-500"

# If the first argument is not empty, use it as the array jobs range
# Otherwise, use the default range
array_jobs="10201,10340,10571,1058,1174,12279,12828,13489,1443,1573,1641,16662,1708,17293,1756,18568,1952,20001,20027,21332,22126,22320,22797,22837,23365,240,25072,2810,2885,3039,3210,322,3264,327,3406,3416,3428,3816,3980,4273,4451,46,597,6103,6245,6716,8453,9620"

echo sbatch --job-name "PHYLO_IQTREE_LOCI" --output "${root}/output_data/slurm_logs/04_Phylogeny_And_Admixture/phylo_run_iqtree_loci_%a.log" -a "${array_jobs}" --export=root=${root} "${root}/code/04_Phylogeny_and_Admixture/run_iqtree_gene_loci_ind_longtime.sb"
sbatch --job-name "PHYLO_IQTREE_LOCI" --output "${root}/output_data/slurm_logs/04_Phylogeny_And_Admixture/phylo_run_iqtree_loci_%a.log" -a "${array_jobs}" --export=root=${root} "${root}/code/04_Phylogeny_and_Admixture/run_iqtree_gene_loci_ind_longtime.sb"
