#!/bin/bash

## 02_whatshap_submit.sh
## Feb 2023
## JRG
# This script creates runs whatshap on sample batches located in the metadata folder
# e.g. bash 03_whatshap_submit.sh sample_names_b_0.txt

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

outdir=${root}/input_data/08_Phasing/
mkdir -p ${outdir}

mkdir -p ${root}/output_data/slurm_logs/08_Phasing/

while read SAMPLE
do
  echo "$SAMPLE"
  mkdir -p ${outdir}/$SAMPLE

  echo sbatch --job-name ${SAMPLE}_WHATSHAP --output ${root}/output_data/slurm_logs/08_Phasing/${SAMPLE}_SQUISH_OUTPUTS.slurm.log --export=root=${root} ${root}/code/08_Phasing/squish_outputs.sb $SAMPLE
  sbatch --job-name ${SAMPLE}_WHATSHAP --output ${root}/output_data/slurm_logs/08_Phasing/${SAMPLE}_SQUISH_OUTPUTS.slurm.log --export=root=${root} ${root}/code/08_Phasing/squish_outputs.sb $SAMPLE

done < "$1"