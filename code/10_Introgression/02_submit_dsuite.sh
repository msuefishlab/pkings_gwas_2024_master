## 02_submit_dsuite.sh
## August 2025
## JRG
## Submits to queue to Run Dsuite for whole genome.

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

mkdir -p ${root}/output_data/slurm_logs/10_TWISST/

outdir=${root}/output_data/10_TWISST/
indir=${root}/input_data/10_TWISST/

mkdir -p ${outdir}

for chr in {1..25}; do
	echo sbatch --job-name "DSUITE_chr${chr}" \
		   --output "${root}/output_data/slurm_logs/10_TWISST/dsuite_${chr}.log" \
		   --export=chr=${chr},root=${root} \
		   "${root}/code/10_TWISST/run_dsuite.sb"

	sbatch --job-name "DSUITE_chr${chr}" \
		   --output "${root}/output_data/slurm_logs/10_TWISST/dsuite_${chr}.log" \
		   --export=chr=${chr},root=${root} \
		   "${root}/code/10_TWISST/run_dsuite.sb"
done