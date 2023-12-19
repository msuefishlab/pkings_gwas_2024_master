## run_upload.sh
## 12-03-23
## JRG
## Submits to queue for variant table construction

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source $root/"pkings_gwas.env"
nfiles=$(wc -l < $1)
inputdata=$1

mkdir -p $root"/output_data/slurm_logs"/03_QC/

echo sbatch --job-name "QC_MAKE_TABLE" --output $root"/output_data/slurm_logs"/03_QC/"qc_make_table.log" --export=root=$root,vcf_file=$vcf_file $root/code/03_QC/get_variant_table.sb
sbatch --job-name "QC_MAKE_TABLE" --output $root"/output_data/slurm_logs"/03_QC/"qc_make_table.log" --export=root=$root,vcf_file=$vcf_file $root/code/03_QC/get_variant_table.sb
