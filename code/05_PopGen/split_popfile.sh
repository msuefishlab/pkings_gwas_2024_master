#!/bin/bash

# Check if the input file and output directory are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_file> <output_directory>"
    exit 1
fi

input_file="$1"
output_dir="$2"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Read the file line by line
while IFS=$'\t' read -r col1 col2 || [ -n "$col1" ]; do
    # Write the first column value to a file in the specified output directory
    echo "$col1" >> "${output_dir}/${col2}.txt"
done < "$input_file"
