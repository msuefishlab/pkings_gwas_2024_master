#!/bin/bash

## Run by bash code/09_RNASeq/02_cat_by_specno.sh input_data/09_RNASeq/rnaseq_metadata.txt

root="$(git rev-parse --show-toplevel)"
source "$root/pkings_gwas.env"

# Define the input file

outdir="${root}/output_data/09_RNASeq/rmats"  # Define the output director
slurmoutdir=${root}/output_data/slurm_outputs/09_RNASeq

mkdir -p $outdir/gw1
mkdir -p $outdir/gw2

mkdir -p $slurmoutdir


echo sbatch --output=${slurmoutdir}/RMATS_post_GW1.out --export=BG1="${root}/input_data/09_RNASeq/g1_bp.txt",BG2="${root}/input_data/09_RNASeq/g1_tp.txt",OUTDIR=${outdir}/gw1, ${root}/code/09_RNASeq/run_rmats.sb;
sbatch --output=${slurmoutdir}/RMATS_post_GW1.out --export=BG1="${root}/input_data/09_RNASeq/g1_bp.txt",BG2="${root}/input_data/09_RNASeq/g1_tp.txt",OUTDIR=${outdir}/gw1, ${root}/code/09_RNASeq/run_rmats.sb;

echo sbatch --output=${slurmoutdir}/RMATS_post_GW2.out --export=BG1="${root}/input_data/09_RNASeq/g2_bp.txt",BG2="${root}/input_data/09_RNASeq/g2_tp.txt",OUTDIR=${outdir}/gw2, ${root}/code/09_RNASeq/run_rmats.sb;
sbatch --output=${slurmoutdir}/RMATS_post_GW2.out --export=BG1="${root}/input_data/09_RNASeq/g2_bp.txt",BG2="${root}/input_data/09_RNASeq/g2_tp.txt",OUTDIR=${outdir}/gw2, ${root}/code/09_RNASeq/run_rmats.sb;