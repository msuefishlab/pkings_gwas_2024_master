## 01_Run_SNP_Phylo.sh
## Dec 2023
## JRG
## Submits to queue for Phylogeny Reconstruction

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

mkdir -p ${root}/output_data/slurm_logs/04_Phylogeny/

echo sbatch --job-name "PHYLO_FILTER_SNPS" --output ${root}"/output_data/slurm_logs/04_Phylogeny/phylo_filter_snps.log" --export=root=${root} ${root}/code/04_Phylogeny/filter_snps.sb
sbatch --job-name "PHYLO_FILTER_SNPS" --output ${root}"/output_data/slurm_logs/04_Phylogeny/phylo_filter_snps.log" --export=root=${root} ${root}/code/04_Phylogeny/filter_snps.sb
