#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

indir=${root}/input_data/04_Phylogeny_And_Admixture/

# Loop over each file in the directory
for FILE in "${indir}"/*.args; do
	bn=$(basename "$FILE")
	echo sbatch --job-name "PREP_POMO" --output ${root}"/output_data/slurm_logs/04_Phylogeny_And_Admixture/phylo_prep_pomo_${bn}.log" --export=root=${root},pop=${bn} ${root}/code/04_Phylogeny_and_Admixture/prep_for_pomo.sb
	sbatch --job-name "PREP_POMO" --output ${root}"/output_data/slurm_logs/04_Phylogeny_And_Admixture/phylo_prep_pomo_${bn}.log" --export=root=${root},pop=${bn} ${root}/code/04_Phylogeny_and_Admixture/prep_for_pomo.sb

done