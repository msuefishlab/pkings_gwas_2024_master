## run_merge.sh
## 12-03-23
## JRG
## Merges interval split tab delimited files for QC Analysis

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

mkdir -p $root"/output_data/slurm_logs"/03_QC/

echo sbatch --job-name "GET_GQs" --output ${root}"/output_data/slurm_logs"/03_QC/"qc_get_gqs.log" --export=root=${root} ${root}/code/03_VCF_QC_Filter/get_GQs.sb
sbatch --job-name "GET_GQs" --output ${root}"/output_data/slurm_logs"/03_QC/"qc_get_gqs.log" --export=root=${root} ${root}/code/03_VCF_QC_Filter/get_GQs.sb
