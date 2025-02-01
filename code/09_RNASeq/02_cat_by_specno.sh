#!/bin/bash

## Run by bash code/09_RNASeq/02_cat_by_specno.sh input_data/09_RNASeq/rnaseq_metadata.txt

root="$(git rev-parse --show-toplevel)"
source "$root/pkings_gwas.env"

# Define the input file
input_file=$1

outdir="${root}/output_data/09_RNASeq/raw_reads"  # Define the output directory

mkdir -p "${outdir}"

# Create associative arrays to store files to concatenate for each SpecNo
declare -A r1_files
declare -A r2_files

# Skip the header and read the file line by line using a file descriptor
{
  read -r header # Skip the header
  while IFS=$'\t' read -r specno population r1 r2 phenotype; do
    # Append the R1 and R2 file paths to the respective associative arrays
    r1_files["$specno"]+="${root}/input_data/09_RNASeq/reads/$r1 "
    r2_files["$specno"]+="${root}/input_data/09_RNASeq/reads/$r2 "
  done
} < "$input_file"

# Iterate over each SpecNo and concatenate the files
for specno in "${!r1_files[@]}"; do
  echo "Concatenating files for SpecNo: $specno"
  
  # Ensure there are files to concatenate
  if [[ -n ${r1_files["$specno"]} ]]; then
    #echo cat ${root}/input_data/09_RNASeq/reads/${r1_files["$specno"]}
    cat ${r1_files["$specno"]} > "${outdir}/${specno}_R1_combined.fastq.gz"
  fi

  if [[ -n ${r2_files["$specno"]} ]]; then
    #echo cat ${root}/input_data/09_RNASeq/reads/${r2_files["$specno"]}
    cat ${r2_files["$specno"]} > "${outdir}/${specno}_R2_combined.fastq.gz"
  fi
done

echo "Concatenation complete."
