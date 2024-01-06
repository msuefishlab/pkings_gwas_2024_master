#!/bin/bash

## 02_run_rel_mat.sh
## January 2024
## JRG
# This script submits a set of individuals specified in arg1 
# e.g.  bash code/06_Association/02_run_rel_mat.sh APA_BAM_ONLY_BP1_TP0_WOBBLE9
# note that output is names based on the keep file.

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

mkdir -p ${root}/output_data/slurm_logs/06_Assocation/

OUTNAME=$1

echo sbatch --job-name $OUTNAME"_MAKE_RELMAT" --output ${root}"/output_data/slurm_logs/06_Assocation/ASSOC_"${OUTNAME}"_MAKE_RELMAT.slurm.log" --export=root=${root},OUTNAME=${OUTNAME} ${root}/code/06_Association/relate_matrix.sb
sbatch --job-name $OUTNAME"_MAKE_RELMAT" --output ${root}"/output_data/slurm_logs/06_Assocation/ASSOC_"${OUTNAME}"_MAKE_RELMAT.slurm.log" --export=root=${root},OUTNAME=${OUTNAME} ${root}/code/06_Association/relate_matrix.sb
