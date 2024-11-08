## copy_files.sh
## 11-08-24
## JRG
## Downloads RNASeq Files from AWS

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source $root/"pkings_gwas.env"
nfiles=$(wc -l < $1)
inputdata=$1

mkdir -p $root"/output_data/slurm_logs"/09_RNASeq/

echo sbatch --job-name "copy_RNAseq" --output $root"/output_data/slurm_logs"/09_RNASeq/"copy_rnaseq_data_%a.log" -a 0-$nfiles --export=root=$root,inputdata=$inputdata $root/code/09_RNASeq/submit_copy.sb
sbatch --job-name "copy_RNAseq" --output $root"/output_data/slurm_logs"/09_RNASeq/"copy_rnaseq_data_%a.log" -a 0-$nfiles --export=root=$root,inputdata=$inputdata $root/code/09_RNASeq/submit_copy.sb
