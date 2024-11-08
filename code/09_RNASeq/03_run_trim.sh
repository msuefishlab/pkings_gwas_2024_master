#!/bin/bash

## Run by bash code/09_RNASeq/02_cat_by_specno.sh input_data/09_RNASeq/rnaseq_metadata.txt

root="$(git rev-parse --show-toplevel)"
source "$root/pkings_gwas.env"

# Define the input file
input_file=$1

outdir="${root}/output_data/09_RNASeq/trimmed_reads"  # Define the output director


# Skip the header and read the file line by line using a file descriptor
{
  read -r header # Skip the header
  while IFS=$'\t' read -r specno population r1 r2 phenotype; do
    mkdir -p ${outdir}/${specno};
    echo sbatch --output=${sample}-trim-%j.out --export=READ1="${i}",READ2="${i//R1/R2}",OUTNAME=$SCRATCH/pkings_trimmed/${sample},SAMPLE=${sample} run_trim.sb;
    sbatch --output=${sample}-trim-%j.out --export=READ1="${i}",READ2="${i//R1/R2}",OUTNAME=$SCRATCH/pkings_trimmed/${sample},SAMPLE=${sample} run_trim.sb;
  done
} < "$input_file"
