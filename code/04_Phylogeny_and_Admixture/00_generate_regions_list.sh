#!/bin/bash

# Check if a file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <GFF_file>"
    exit 1
fi

# Check if the file exists
if [ ! -f "$1" ]; then
    echo "File not found: $1"
    exit 1
fi

# Run the awk command to extract and format gene information
awk '$3 == "gene" { 
  split($9, a, ";"); 
  for (i in a) {
    if (a[i] ~ /^ID=/) { 
      split(a[i], id, "="); 
      # Remove "gene-" prefix if present
      geneID=gensub(/^gene-/, "", "g", id[2]); 
    }
  }
  printf "%s@%s:%s-%s\n", geneID, $1, $4, $5 
}' "$1"
