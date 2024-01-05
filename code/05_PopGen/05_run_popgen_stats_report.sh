## 05_run_popgen_stats_report.sh
## 01-20204
## JRG
## Submits jobs to slurm for constructing Overall Popgen Statistics Report Files

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

mkdir -p $root"/output_data/slurm_logs"/05_PopGen/

echo sbatch --job-name "RUN_POPGEN_REPORT" --output ${root}"/output_data/slurm_logs"/05_PopGen/"POPGEN_report.log" --export=root=${root} ${root}/code/05_PopGen/render_popgen_report.sb
