#!/bin/bash

root="$(git rev-parse --show-toplevel)"

source ${root}/"pkings_gwas.env"

CMD="cd ${root}; module purge; module load SciPy-bundle/2023.11-gfbf-2023b; python3 code/03_VCF_QC_Filter/13_generate_missing_table.py"

echo sbatch --time=04:00:00 -c 8 -J GENERATE_MISSING_TABLE -o ${root}"/output_data/slurm_logs/03_QC/QC_generate_missing.log" --export=root=${root} --mem 12000 --wrap="$CMD"
sbatch --time=04:00:00 -c 8 -J GENERATE_MISSING_TABLE -o ${root}"/output_data/slurm_logs/03_QC/QC_generate_missing.log" --export=root=${root} --mem 12000 --wrap="$CMD"
