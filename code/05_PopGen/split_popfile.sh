#!/bin/bash

# Check if the input file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <filename>"
    exit 1
fi

# Read the file line by line
while IFS=$'\t' read -r col1 col2 || [ -n "$col1" ]; do
    # Write the first column value to a file named by the second column value
    echo "$col1" >> "${col2}.txt"
done < "$1"
