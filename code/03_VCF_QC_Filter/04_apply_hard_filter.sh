## run_upload.sh
## 12-03-23
## JRG
## Submits to queue for variant table construction

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

interval_files=($(find "${root}/input_data/01_Terra/chrom-interval-files/" "${root}/input_data/01_Terra/remaining-interval-files/" -type f))

nintervals=${#interval_files[@]}

mkdir -p $root"/output_data/slurm_logs"/03_QC/

echo sbatch --job-name "QC_FILTER_SNPS" --output ${root}"/output_data/slurm_logs"/03_QC/"qc_filter_snps_%a.log" -a 0-$nintervals --export=root=${root},vartype="SNP" ${root}/code/03_VCF_QC_Filter/apply_hard_filter.sb
sbatch --job-name "QC_FILTER_SNPS" --output ${root}"/output_data/slurm_logs"/03_QC/"qc_filter_snps_%a.log" -a 0-$nintervals --export=root=${root},vartype="SNP" ${root}/code/03_VCF_QC_Filter/apply_hard_filter.sb

echo sbatch --job-name "QC_FILTER_INDELS" --output ${root}"/output_data/slurm_logs"/03_QC/"qc_filter_indels_%a.log" -a 0-$nintervals --export=root=${root},vartype="INDEL" ${root}/code/03_VCF_QC_Filter/apply_hard_filter.sb
sbatch --job-name "QC_FILTER_SNPS" --output ${root}"/output_data/slurm_logs"/03_QC/"qc_filter_indels_%a.log" -a 0-$nintervals --export=root=${root},vartype="INDEL" ${root}/code/03_VCF_QC_Filter/apply_hard_filter.sb