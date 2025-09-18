#!/bin/bash

## Run by bash code/09_RNASeq/02_cat_by_specno_tissue.sh input_data/09_RNASeq/rnaseq_metadata.txt

root="$(git rev-parse --show-toplevel)"
source "$root/pkings_gwas.env"

# Define the input file
input_file=$1

readsdir="${root}/output_data/09_RNASeq/trimmed_reads_somatic"  # Define the directory for trimmed reads
outdir="${root}/output_data/09_RNASeq/aligned_reads"  # Define the directory for aligned reads
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
      mkdir -p "${outdir}/${key}"
      echo sbatch --output="${slurmoutdir}/align_${key}.out" --export=READ1="${readsdir}/${key}/tr_${key}_1P.fastq.gz",READ2="${readsdir}/${key}/tr_${key}_2P.fastq.gz",OUTDIR="${outdir}/${key}",SAMPLE="${key}" "${root}/code/09_RNASeq/run_star_rnaseq.sb"
      sbatch --output="${slurmoutdir}/align_${key}.out" --export=READ1="${readsdir}/${key}/tr_${key}_1P.fastq.gz",READ2="${readsdir}/${key}/tr_${key}_2P.fastq.gz",OUTDIR="${outdir}/${key}",SAMPLE="${key}" "${root}/code/09_RNASeq/run_star_rnaseq.sb"
    fi
  done
} < "$input_file"
