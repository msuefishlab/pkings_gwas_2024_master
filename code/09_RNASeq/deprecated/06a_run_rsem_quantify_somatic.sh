#!/bin/bash

## Run by bash code/09_RNASeq/02_cat_by_specno_tissue.sh input_data/09_RNASeq/rnaseq_metadata.txt

root="$(git rev-parse --show-toplevel)"
source "$root/pkings_gwas.env"

# Define the input file
input_file=$1

indir="${root}/output_data/09_RNASeq/aligned_reads_somatic"  # Directory for aligned reads
outdir="${root}/output_data/09_RNASeq/quantification"  # Directory for quantification output
slurmoutdir="${root}/output_data/slurm_outputs/09_RNASeq"

mkdir -p "$outdir"
mkdir -p "$slurmoutdir"

# Skip the header and read the file line by line using a file descriptor
{
  read -r header # Skip the header
  declare -A seen_keys # Declare an associative array to track unique specno_tissue keys

  while IFS=$'\t' read -r specno tissue r1 r2; do
    key="${specno}_${tissue}" # Create the unique key based on SpecNo and Tissue
    
    if [[ -z "${seen_keys[$key]}" ]]; then
      # If key is not yet seen
      seen_keys[$key]=1 # Mark key as seen
      echo sbatch --output="${slurmoutdir}/RSEM_${key}.out" --export=INDIR="${indir}/${key}",OUTDIR="${outdir}/${key}",SAMPLE="${key}" "${root}/code/09_RNASeq/run_rsem.sb"
      sbatch --output="${slurmoutdir}/RSEM_${key}.out" --export=INDIR="${indir}/${key}",OUTDIR="${outdir}/${key}",SAMPLE="${key}" "${root}/code/09_RNASeq/run_rsem.sb"
    fi
  done
} < "$input_file"
