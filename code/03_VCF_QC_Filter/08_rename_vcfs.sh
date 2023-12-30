#!/bin/bash

## 08_rename_vcfs.sh
## Dec 2023
## JRG
# This script renames assembly scaffolds based on the chr_name_map file specified in the env file

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"
scratchdir=${scratch_store}/output_data/03_QC/
outdir=${root}/output_data/03_QC/


echo sbatch --job-name "RENAME_CHRS" --output ${root}"/output_data/slurm_logs"/03_QC/"qc_rename_%a.log " -a 0-2 --export=root=${root} ${root}/code/03_VCF_QC_Filter/rename_vcfs.sb
sbatch --job-name "RENAME_CHRS" --output ${root}"/output_data/slurm_logs"/03_QC/"qc_rename_%a.log " -a 0-2 --export=root=${root} ${root}/code/03_VCF_QC_Filter/rename_vcfs.sb
