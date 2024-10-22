#!/bin/bash

## 02_whatshap_submit.sh
## Feb 2023
## JRG
# This script creates runs whatshap on sample batches located in the metadata folder
# e.g. bash 03_whatshap_submit.sh sample_names_b_0.txt

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

outdir=${root}/input_data/07_Phasing/
mkdir -p ${outdir}

mkdir -p ${root}/output_data/slurm_logs/07_Phasing/

while read SAMPLE
do
  echo "$SAMPLE"
  mkdir -p ${outdir}/$SAMPLE

  echo sbatch --job-name ${SAMPLE}_WHATSHAP --array=0-25 --output ${root}/output_data/slurm_logs/07_Phasing/${SAMPLE}_WHATSHAP.slurm_%a.log --export=root=${root} ${root}/code/07_Phasing/submit_whatshap.sb $SAMPLE
  sbatch --job-name ${SAMPLE}_WHATSHAP --array=0-25 --output ${root}/output_data/slurm_logs/07_Phasing/${SAMPLE}_WHATSHAP.slurm_%a.log --export=root=${root} ${root}/code/07_Phasing/submit_whatshap.sb $SAMPLE

done < "$1"