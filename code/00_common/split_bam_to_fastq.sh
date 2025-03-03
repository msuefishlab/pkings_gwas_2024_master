#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source "$root/pkings_gwas.env"

# Define the input file

outdir="${root}/output_data/11_NCBI_Upload" # Define the output director
slurmoutdir=${root}/output_data/slurm_logs/11_NCBI_Upload

mkdir -p $outdir

mkdir -p $slurmoutdir



for bam in /mnt/research/efish/incoming_new_assemblies/MIC4273/NewHap2/delivery/bams/*.bam; do
	bn_bam=$(basename $bam)
    echo sbatch --output=${slurmoutdir}/${bn_bam}_convert_to_fq.out --export=OUTDIR=${outdir},BAM=${bam}, ${root}/code/00_common/split_bam.sb;
	sbatch --output=${slurmoutdir}/${bn_bam}_convert_to_fq.out --export=OUTDIR=${outdir},BAM=${bam}, ${root}/code/00_common/split_bam.sb;
done