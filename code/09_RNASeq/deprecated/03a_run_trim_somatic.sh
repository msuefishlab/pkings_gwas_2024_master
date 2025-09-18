#!/bin/bash

## Run by bash code/09_RNASeq/02_cat_by_specno_tissue.sh input_data/09_RNASeq/rnaseq_metadata.txt

root="$(git rev-parse --show-toplevel)"
source "$root/pkings_gwas.env"

# Define the input file
input_file=$1

readsdir="${root}/output_data/09_RNASeq/raw_reads"  # Define the directory for raw reads
outdir="${root}/output_data/09_RNASeq/trimmed_reads"  # Define the directory for trimmed reads
slurmoutdir="${root}/output_data/slurm_outputs/09_RNASeq"

mkdir -p "$outdir"
mkdir -p "$slurmoutdir"

{
  read -r header # Skip the header
  declare -A seen_keys # Declare an associative array to track unique specno_tissue keys

  while IFS=$'\t' read -r specno tissue r1 r2; do
    key="${specno}_${tissue}" # Create the unique key based on SpecNo and Tissue
    
    if [[ -z "${seen_keys[$key]}" ]]; then
      # If key is not yet seen
      seen_keys[$key]=1 # Mark key as seen
      mkdir -p "${outdir}/${key}"
      echo sbatch --output="${slurmoutdir}/trim_${key}.out" --export=READ1="${readsdir}/${key}_R1_combined.fastq.gz",READ2="${readsdir}/${key}_R2_combined.fastq.gz",OUTNAME="${outdir}/${key}",SAMPLE="${key}" "${root}/code/09_RNASeq/run_trim.sb"
      sbatch --output="${slurmoutdir}/trim_${key}.out" --export=READ1="${readsdir}/${key}_R1_combined.fastq.gz",READ2="${readsdir}/${key}_R2_combined.fastq.gz",OUTNAME="${outdir}/${key}",SAMPLE="${key}" "${root}/code/09_RNASeq/run_trim.sb"
    fi
  done
} < "$input_file"
