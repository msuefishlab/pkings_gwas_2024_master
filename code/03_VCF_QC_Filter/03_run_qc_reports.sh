## 03_run_qc_reports.sh
## 12-03-23
## JRG
## Submits jobs to slurm for constructing QC Report Files

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

mkdir -p $root"/output_data/slurm_logs"/03_QC/

echo sbatch --job-name "RUN_QC_REPORTS" --output ${root}"/output_data/slurm_logs"/03_QC/"qc_report.log" --export=root=${root} ${root}/code/03_VCF_QC_Filter/Render_QC_Report.sb
sbatch --job-name "RUN_QC_REPORTS" --output ${root}"/output_data/slurm_logs"/03_QC/"qc_report.log" --export=root=${root} ${root}/code/03_VCF_QC_Filter/Render_QC_Report.sb
