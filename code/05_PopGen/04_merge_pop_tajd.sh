## 03_merge_tajd.sh
## 12-03-23
## JRG
## Merges chromsome split Tajima's D Files into Population TajD Files

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

# Define the output file
output_file=${root}/output_data/05_PopGen/${popgen_prefix}.tajd.stats.final.txt

# Check if the output file already exists and remove it to start fresh
if [ -f "$output_file" ]; then
    rm "$output_file"
fi

# Header flag
header_written=false

# Loop through each .stats.txt file
for file in ${root}/output_data/05_PopGen/*.tajd.stats.txt; do
    # Extract the population name from the filename
    pop_name=$(basename "$file" .tajd.stats.txt)

    # Check if header is written to output file
    if [ "$header_written" = false ] ; then
        # Write the header with the new 'pop' column
        echo -e "pop1\tCHROM\tBIN_START\tN_SNPS\tTajimaD" >> "$output_file"
        header_written=true
    fi

    # Add the population name as the first column, skip the second line (blank line)
    awk -v pop="$pop_name" 'NR>2 {print pop, $0}' OFS="\t" "$file" >> "$output_file"
done

echo "Combined file created: $output_file"
