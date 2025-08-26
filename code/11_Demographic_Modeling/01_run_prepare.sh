#!/bin/bash

## 01_run_prep_data.sh
## January 2024
## JRG
# This script prepares a VCF files for a set of individuals specified in arg1 with phenotypes in arg2 for relationship matrix construction and association analysis.
# e.g. bash code/06_Association/01_run_prep_data.sh APA_BAM_ONLY_BP0_TP1_WOB9

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

mkdir -p ${root}/output_data/slurm_logs/11_Demographic_Modeling/

echo sbatch --job-name $OUTNAME"_PREP_DATA" --output ${root}"/output_data/slurm_logs/11_Demographic_Modeling/prep_data.log" --export=root=${root} ${root}/code/11_Demographic_Modeling/prepare.sb
sbatch --job-name $OUTNAME"_PREP_DATA" --output ${root}"/output_data/slurm_logs/11_Demographic_Modeling/prep_data.log" --export=root=${root} ${root}/code/11_Demographic_Modeling/prepare.sb
