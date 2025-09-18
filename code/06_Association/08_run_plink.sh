#!/bin/bash

## 03_run_gemma.sh
## January 2024
## JRG
# This script submits a set of individuals specified in arg1 to run a gemma association test
# e.g.  bash code/06_Association/03_run_gemma.sh APA_BAM_ONLY_BP1_TP0_WOBBLE9
# note that output is names based on the keep file.

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

mkdir -p ${root}/output_data/slurm_logs/06_Assocation/

OUTNAME=$1

echo sbatch --job-name $OUTNAME"_PLINK" --output ${root}"/output_data/slurm_logs/06_Assocation/ASSOC_"${OUTNAME}"_RUN_PLINK.log" --export=root=${root},OUTNAME=${OUTNAME} ${root}/code/06_Association/plink.sb
sbatch --job-name $OUTNAME"_PLINK" --output ${root}"/output_data/slurm_logs/06_Assocation/ASSOC_"${OUTNAME}"_RUN_PLINK.slurm.log" --export=root=${root},OUTNAME=${OUTNAME} ${root}/code/06_Association/plink.sb