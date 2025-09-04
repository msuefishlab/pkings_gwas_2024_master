#!/bin/bash

## 03_run_gemma.sh
## January 2024
## JRG
# This script submits a set of individuals specified in arg1 to run a gemma association test
# e.g.  bash code/06_Association/03_run_gemma.sh APA_BAM_ONLY_BP1_TP0_WOBBLE9
# note that output is names based on the keep file.

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

mkdir -p ${root}/output_data/slurm_logs/14_Wobble_Analysis/

OUTNAME=$1

echo sbatch --job-name $OUTNAME"_GEMMA_COVARS" --output ${root}"/output_data/slurm_logs/14_Wobble_Analysis/ASSOC_"${OUTNAME}"_RUN_GEMMA_COVARS.log" --export=root=${root},OUTNAME=${OUTNAME} ${root}/code/14_Wobble_Analysis/gemma_w_covars.sb
sbatch --job-name $OUTNAME"_GEMMA_COVARS" --output ${root}"/output_data/slurm_logs/14_Wobble_Analysis/ASSOC_"${OUTNAME}"_RUN_GEMMA_COVARS.slurm.log" --export=root=${root},OUTNAME=${OUTNAME} ${root}/code/14_Wobble_Analysis/gemma_w_covars.sb