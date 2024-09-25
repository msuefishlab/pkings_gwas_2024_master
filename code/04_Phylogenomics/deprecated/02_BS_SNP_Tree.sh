## 02_BS_SNP_Tree.sh
## Jan 2024
## JRG
## Submits to queue for Phylogeny Bootstrap

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

mkdir -p ${root}/output_data/slurm_logs/04_Phylogeny/

echo sbatch --job-name "PHYLO_BS_PHYLO" --output ${root}"/output_data/slurm_logs/04_Phylogeny/phylo_bootstrap.log" --export=root=${root} ${root}/code/04_Phylogeny/bootstrap_tree.sb
sbatch --job-name "PHYLO_BS_PHYLO" --output ${root}"/output_data/slurm_logs/04_Phylogeny/phylo_bootstrap.log" --export=root=${root} ${root}/code/04_Phylogeny/bootstrap_tree.sb
