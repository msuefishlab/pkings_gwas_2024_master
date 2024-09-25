#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

indir=${root}/input_data/04_Phylogeny_And_Admixture/


echo sbatch --job-name "RUN_POMO" --output ${root}"/output_data/slurm_logs/04_Phylogeny_And_Admixture/phylo_run_pomo.log" --export=root=${root} ${root}/code/04_Phylogeny_and_Admixture/run_iqtree_pomo.sb
sbatch --job-name "RUN_POMO" --output ${root}"/output_data/slurm_logs/04_Phylogeny_And_Admixture/phylo_run_pomo.log" --export=root=${root} ${root}/code/04_Phylogeny_and_Admixture/run_iqtree_pomo.sb

