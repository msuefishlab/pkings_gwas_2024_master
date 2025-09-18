#!/bin/bash

## Run by bash code/09_RNASeq/02_cat_by_specno_tissue.sh input_data/09_RNASeq/rnaseq_metadata.txt

root="$(git rev-parse --show-toplevel)"
source "$root/pkings_gwas.env"

# Define the input file
input_file=$1

outdir="${root}/output_data/09_RNASeq/raw_reads"  # Define the output directory

mkdir -p "${outdir}"

# Create associative arrays to store files to concatenate for each SpecNo + Tissue
declare -A r1_files
declare -A r2_files

# Skip the header and read the file line by line using a file descriptor
{
  read -r header # Skip the header
  while IFS=$'\t' read -r specno tissue r1 r2; do
    key="${specno}_${tissue}"  # Use SpecNo and Tissue as the key
    # Append the R1 and R2 file paths to the respective associative arrays
    r1_files["$key"]+="${root}/input_data/09_RNASeq/reads/$r1 "
    r2_files["$key"]+="${root}/input_data/09_RNASeq/reads/$r2 "
  done
} < "$input_file"

# Iterate over each SpecNo_Tissue combination and concatenate the files
for key in "${!r1_files[@]}"; do
  echo "Concatenating files for SpecNo_Tissue: $key"
  
  # Ensure there are files to concatenate
  if [[ -n ${r1_files["$key"]} ]]; then
    cat ${r1_files["$key"]} > "${outdir}/${key}_R1_combined.fastq.gz"
  fi

  if [[ -n ${r2_files["$key"]} ]]; then
    cat ${r2_files["$key"]} > "${outdir}/${key}_R2_combined.fastq.gz"
  fi
done

echo "Concatenation complete."
