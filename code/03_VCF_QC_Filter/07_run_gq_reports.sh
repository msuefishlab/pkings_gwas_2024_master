## run_merge.sh
## 12-03-23
## JRG
## Submits the GQ table for constructing a GQ Report File

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

mkdir -p $root"/output_data/slurm_logs"/03_QC/

echo sbatch --job-name "RENDER_GQ_REPORTS" --output ${root}"/output_data/slurm_logs"/03_QC/"qc_gq_report.log" --export=root=${root} ${root}/code/03_VCF_QC_Filter/Render_GQ_Report.sb
sbatch --job-name "RENDER_GQ_REPORTS" --output ${root}"/output_data/slurm_logs"/03_QC/"qc_gq_report.log" --export=root=${root} ${root}/code/03_VCF_QC_Filter/Render_GQ_Report.sb
