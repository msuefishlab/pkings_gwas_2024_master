## Submit_TWISST.sh
## July 2025
## JRG
## Submits to queue to Run TWISST for chromosomes containing peaks.

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

mkdir -p ${root}/output_data/slurm_logs/10_TWISST/

outdir=${root}/output_data/10_TWISST/
indir=${root}/input_data/10_TWISST/

mkdir -p ${outdir}

# Read chromosome list into an array
chr_list=()
while IFS= read -r line; do
    chr_list+=("$line")
done < "${indir}/chr_list.txt"

for chr in "${chr_list[@]}"; do
	echo sbatch --job-name "TWISST_${chr}" \
		   --output "${root}/output_data/slurm_logs/10_TWISST/prep_twisst_${chr}.log" \
		   --export=chr=${chr},root=${root} \
		   "${root}/code/10_TWISST/prep_twisst.sb"

	sbatch --job-name "TWISST_${chr}" \
		   --output "${root}/output_data/slurm_logs/10_TWISST/prep_twisst_${chr}.log" \
		   --export=chr=${chr},root=${root} \
		   "${root}/code/10_TWISST/prep_twisst.sb"
done