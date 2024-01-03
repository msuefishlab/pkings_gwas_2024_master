## 01_submit_variant_table.sh
## 12-03-23
## JRG
## Submits to queue for variant table construction

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

interval_files=($(find "${root}/input_data/01_Terra/chrom-interval-files/" "${root}/input_data/01_Terra/remaining-interval-files/" -type f))

nintervals=${#interval_files[@]}

mkdir -p $root"/output_data/slurm_logs"/05_PopGen/

echo sbatch --job-name "POPGEN_TAJD" --output ${root}"/output_data/slurm_logs"/05_PopGen/"POPGEN_TAJD_%a.log" -a 0-$nintervals --export=root=${root} ${root}/code/05_PopGen/tajd.sb
sbatch --job-name "POPGEN_TAJD" --output ${root}"/output_data/slurm_logs"/05_PopGen/"POPGEN_TAJD_%a.log" -a 0-$nintervals --export=root=${root} ${root}/code/05_PopGen/tajd.sb