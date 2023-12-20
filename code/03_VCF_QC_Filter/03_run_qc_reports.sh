## run_merge.sh
## 12-03-23
## JRG
## Merges interval split tab delimited files for QC Analysis

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

mkdir -p $root"/output_data/slurm_logs"/03_QC/

echo sbatch --job-name "RUN_QC_REPORTS" --output ${root}"/output_data/slurm_logs"/03_QC/"qc_report.log" --export=root=${root} ${root}/code/03_VCF_QC_Filter/Render_QC_Report.sb
sbatch --job-name "RUN_QC_REPORTS" --output ${root}"/output_data/slurm_logs"/03_QC/"qc_report.log" --export=root=${root} ${root}/code/03_VCF_QC_Filter/Render_QC_Report.sb
