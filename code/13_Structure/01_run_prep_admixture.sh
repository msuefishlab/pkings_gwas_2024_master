#!/bin/bash

## 06a_run_prep_admixture.sh
## Jan 2024
## JRG
# This script runs various vcftools to prepare a PLINK file for ADMIXTURE

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

mkdir -p ${root}/output_data/slurm_logs/13_Structure/

echo sbatch --job-name "POPGEN_PREP_ADMIX" --output ${root}"/output_data/slurm_logs/13_Structure/popgen_prep_admixture.log" --export=root=${root} ${root}/code/13_Structure/prep_admixture.sb
sbatch --job-name "POPGEN_PREP_ADMIX" --output ${root}"/output_data/slurm_logs/13_Structure/popgen_prep_admixture.log" --export=root=${root} ${root}/code/13_Structure/prep_admixture.sb
