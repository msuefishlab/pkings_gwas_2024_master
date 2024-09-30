## 01_submit_variant_table.sh
## 12-03-23
## JRG
## Submits to queue for variant table construction

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

intervals=()

# Use a while loop to read the file line by line
while IFS= read -r line; do
    first_column=$(echo "$line" | awk '{print $1}')
    # Append the first column to the array
    intervals+=("$first_column")
done < "${root}/input_data/00_Reference_Genome/${reference}.fai"

nintervals=${#intervals[@]}

mkdir -p $root"/output_data/slurm_logs"/05_PopGen/

echo sbatch --job-name "POPGEN_TAJD" --output ${root}"/output_data/slurm_logs"/05_PopGen/"POPGEN_TAJD_%a.log" -a 0-$nintervals --export=root=${root} ${root}/code/05_PopGen/tajd.sb
sbatch --job-name "POPGEN_TAJD" --output ${root}"/output_data/slurm_logs"/05_PopGen/"POPGEN_TAJD_%a.log" -a 0-$nintervals --export=root=${root} ${root}/code/05_PopGen/tajd.sb