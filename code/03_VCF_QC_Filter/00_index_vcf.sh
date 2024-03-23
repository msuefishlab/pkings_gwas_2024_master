## 00_index_vcf.sh
## 03-23-24
## JRG
## Creates Index for VCF File

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

mkdir -p ${root}/output_data/slurm_logs/03_QC/

echo sbatch --job-name "INDEX_VCF" --output ${root}/output_data/slurm_logs/03_QC/INDEX_VCF.log --export=root=${root} ${root}/code/03_VCF_QC_Filter/index_vcf.sb
sbatch --job-name "INDEX_VCF" --output ${root}/output_data/slurm_logs/03_QC/INDEX_VCF.log --export=root=${root} ${root}/code/03_VCF_QC_Filter/index_vcf.sb
