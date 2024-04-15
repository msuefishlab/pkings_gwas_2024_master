#!/bin/bash

# Directory containing tree files
tree_dir=${root}/

# Directory containing fasta files
fasta_dir="/path/to/fasta_files/"

# Target directory for fasta files
target_dir="/path/to/target_directory/"

# Ensure the target directory exists
mkdir -p "$target_dir"

# Loop through each treefile in the tree file directory
for treefile in "${tree_dir}"*.treefile; do
  # Extract the base name without the extension
  base_name=$(basename "$treefile" .treefile)

  # Construct the corresponding fasta file path
  fasta_file="${fasta_dir}${base_name}.fasta"

  # Check if the fasta file exists and copy it to the target directory
  if [ -f "$fasta_file" ]; then
    cp "$fasta_file" "$target_dir"
  else
    echo "No corresponding fasta file found for $treefile"
  fi
done

