#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source "${root}/pkings_gwas.env"

# Define your array
alignments=("${root}/output_data/04_Phylogeny_And_Admixture/gene_alignments/"*.fasta)

# Get the search file from the first argument
search_file="$1"

# Initialize index
index=-1

# Loop through the array to find the index
for i in "${!alignments[@]}"; do
  if [[ "${alignments[i]}" == *"$search_file" ]]; then
    index=$i
    break
  fi
done

# Check if the file was found
if [[ $index -ge 0 ]]; then
  echo "File found at index: $index"
else
  echo "File not found"
fi
