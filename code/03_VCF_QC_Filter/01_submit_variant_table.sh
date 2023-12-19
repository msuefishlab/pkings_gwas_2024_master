## run_upload.sh
## 12-03-23
## JRG
## Submits to queue for variant table construction

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"


mkdir -p $root"/output_data/slurm_logs"/03_QC/

echo sbatch --job-name "QC_MAKE_TABLE" --output ${root}"/output_data/slurm_logs"/03_QC/"qc_make_table.log" --export=root=${root} ${root}/code/03_QC/get_variant_table.sb
sbatch --job-name "QC_MAKE_TABLE" --output ${root}"/output_data/slurm_logs"/03_QC/"qc_make_table.log" --export=root=${root} ${root}/code/03_QC/get_variant_table.sb
