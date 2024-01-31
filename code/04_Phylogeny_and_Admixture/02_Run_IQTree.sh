## 02_Run_Prep_For_IQTree.sh
## Dec 2023
## JRG
## Submits to queue for Phylogeny Reconstruction

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

mkdir -p ${root}/output_data/slurm_logs/04_Phylogeny_And_Admixture/

echo sbatch --job-name "PHYLO_PREP_IQ_TREE" --output ${root}"/output_data/slurm_logs/04_Phylogeny_And_Admixture/phylo_run_iqtree.log" --export=root=${root} ${root}/code/04_Phylogeny_And_Admixture/run_iqtree.sb
sbatch --job-name "PHYLO_PREP_IQ_TREE" --output ${root}"/output_data/slurm_logs/04_Phylogeny_And_Admixture/phylo_run_iqtree.log" --export=root=${root} ${root}/code/04_Phylogeny_And_Admixture/run_iqtree.sb
