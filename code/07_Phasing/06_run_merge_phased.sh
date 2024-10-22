#!/bin/bash

## 02_whatshap_submit.sh
## Feb 2023
## JRG
# This script runs shapeit on all whatshap vcf files
# e.g. bash 05_run_shapeit.sh

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

outdir=${root}/input_data/07_Phasing/
mkdir -p ${outdir}

mkdir -p ${root}/output_data/slurm_logs/07_Phasing/

echo sbatch --job-name "MERGE_PHASED" --output ${root}/output_data/slurm_logs/07_Phasing/MERGE_PHASED.slurm.log --export=root=${root} ${root}/code/07_Phasing/merge_phased.sb
sbatch --job-name "MERGE_PHASED" --output ${root}/output_data/slurm_logs/07_Phasing/MERGE_PHASED.slurm.log --export=root=${root} ${root}/code/07_Phasing/merge_phased.sb