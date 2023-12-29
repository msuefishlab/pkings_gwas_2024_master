## 06_get_gqs.sh
## 12-03-23
## JRG
## Submits a job to extract GQ values from Biallelic VCF File for QC analysis.

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

mkdir -p $root"/output_data/slurm_logs"/03_QC/

echo sbatch --job-name "GET_GQs" --output ${root}"/output_data/slurm_logs"/03_QC/"qc_get_gqs.log" --export=root=${root} ${root}/code/03_VCF_QC_Filter/get_GQs.sb
sbatch --job-name "GET_GQs" --output ${root}"/output_data/slurm_logs"/03_QC/"qc_get_gqs.log" --export=root=${root} ${root}/code/03_VCF_QC_Filter/get_GQs.sb
