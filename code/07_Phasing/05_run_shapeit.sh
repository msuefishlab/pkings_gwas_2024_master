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

cp ${root}/input_data/01_Terra/chrom_intervals.list ${root}/input_data/07_Phasing/shapeit_regions.txt
echo "final_chunk" >> ${root}/input_data/07_Phasing/shapeit_regions.txt

num_regions=$(cat ${root}/input_data/07_Phasing/shapeit_regions.txt | wc -l)


echo sbatch --job-name ${SAMPLE}_SHAPEIT -a 1-${num_regions} --output ${root}/output_data/slurm_logs/07_Phasing/PHASE_SHAPEIT.slurm_%a.log--export=root=${root} ${root}/code/07_Phasing/run_shapeit.sb
sbatch --job-name ${SAMPLE}_SHAPEIT -a 1-${num_regions} --output ${root}/output_data/slurm_logs/07_Phasing/PHASE_SHAPEIT.slurm_%a.log --export=root=${root} ${root}/code/07_Phasing/run_shapeit.sb