#!/bin/bash

## 01_run_rel_mat.sh
## January 2024
## JRG
# This script submits a set of individuals specified in arg1 with phenotypes in arg2 to make a relationship matrix
# e.g. bash 01_run_rel_mat.sh ./apa_bam_wobs_only.txt ./apa_bam_wobs_only_phenos.txt
# note that output is names based on the keep file.

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

mkdir -p ${root}/output_data/slurm_logs/06_Assocation/

keep_path=$(realpath $1)
pheno_path=$(realpath $2)

keepfilename=$(basename $keep_path)

OUTNAME=${keepfilename%.txt}

echo sbatch --job-name $OUTNAME"_MAKE_RELMAT" --output ${root}"/output_data/slurm_logs/06_Assocation/gemma/ASSOC_"${OUTNAME}"_MAKE_RELMAT.slurm.log" --export=root=${root} ${root}/code/06_Association/relate_matrix.sb $keep_path $pheno_path $OUTNAME
sbatch --job-name $OUTNAME"_MAKE_RELMAT" --output ${root}"/output_data/slurm_logs/06_Assocation/gemma/ASSOC_"${OUTNAME}"_MAKE_RELMAT.slurm.log" --export=root=${root} ${root}/code/06_Association/relate_matrix.sb $keep_path $pheno_path $OUTNAME