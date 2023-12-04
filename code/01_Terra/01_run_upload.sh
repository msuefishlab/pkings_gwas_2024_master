## run_upload.sh
## 12-03-23
## JRG
## Uploads Read Data to Google Cloud/Terra for Processing

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source $root/"pkings_gwas.env"
nfiles=$(wc -l < $1)
inputdata=$1

mkdir -p $root"/output_data/slurm_logs"/01_Terra/

echo sbatch --job-name "copy_data_to_gcloud" --output $root"/output_data/slurm_logs"/01_Terra/"copy_data_to_gcloud_%a.log" -a 0-$nfiles --export=root=$root,inputdata=$inputdata $root/code/01_Terra/submit_copy.sb
sbatch --job-name "copy_data_to_gcloud" --output $root"/output_data/slurm_logs"/01_Terra/"copy_data_to_gcloud_%a.log" -a 0-$nfiles --export=root=$root,inputdata=$inputdata $root/code/01_Terra/submit_copy.sb