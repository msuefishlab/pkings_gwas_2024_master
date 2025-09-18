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

while read -r chr start end peak _; do

echo sbatch --job-name "DSUITE_chr${chr}" \
		   --output "${root}/output_data/slurm_logs/10_TWISST/dsuite_peaks_peak_${peak}.log" \
		   --export=chr=${chr},root=${root},start=${start},end=${end},peak=${peak} \
		   "${root}/code/10_TWISST/run_dsuite_peaks.sb"

	sbatch --job-name "DSUITE_chr${chr}" \
		   --output "${root}/output_data/slurm_logs/10_TWISST/dsuite_peaks_peak_${peak}.log" \
		   --export=chr=${chr},root=${root},start=${start},end=${end},peak=${peak} \
		   "${root}/code/10_TWISST/run_dsuite_peaks.sb"
    
done < ${root}/input_data/10_TWISST/peaks_100kb.bed