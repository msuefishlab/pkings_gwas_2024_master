#!/bin/bash

## 01_run_pixy.sh
## Jan 2024
## JRG
# This script runs pixy to calculate pi, dxy, and fst on VCF file.

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

mkdir -p ${root}/output_data/slurm_logs/05_PopGen/

echo sbatch --job-name "POPGEN_RUN_PIXY" --output ${root}"/output_data/slurm_logs/05_PopGen/popgen_prep_admixture.log" --export=root=${root} ${root}/code/05_PopGen/prep_admixture.sb
sbatch --job-name "POPGEN_RUN_PIXY" --output ${root}"/output_data/slurm_logs/05_PopGen/popgen_prep_admixture.log" --export=root=${root} ${root}/code/05_PopGen/prep_admixture.sb
