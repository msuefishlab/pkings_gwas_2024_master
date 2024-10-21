#!/bin/bash

## 08_rename_vcfs.sh
## Dec 2023
## JRG
# This script renames assembly scaffolds based on the chr_name_map file specified in the env file

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"


echo sbatch --job-name "RENAME_CHRS" --output ${root}/output_data/slurm_logs/08_Phasing/RENAME_PHASED.log --export=root=${root} ${root}/code/08_Phasing/rename_vcfs.sb
sbatch --job-name "RENAME_CHRS" --output ${root}/output_data/slurm_logs/08_Phasing/RENAME_PHASED.log --export=root=${root} ${root}/code/08_Phasing/rename_vcfs.sb
