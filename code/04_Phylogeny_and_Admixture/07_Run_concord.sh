## 02_Run_Prep_For_IQTree.sh
## Dec 2023
## JRG
## Submits to queue for Phylogeny Reconstruction

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

mkdir -p ${root}/output_data/slurm_logs/04_Phylogeny_And_Admixture/

echo sbatch --job-name "PHYLO_ASTRAL" --output ${root}"/output_data/slurm_logs/04_Phylogeny_And_Admixture/phylo_run_concord_extracted.log" --export=root=${root} ${root}/code/04_Phylogeny_and_Admixture/concord.sb
sbatch --job-name "PHYLO_ASTRAL" --output ${root}"/output_data/slurm_logs/04_Phylogeny_And_Admixture/phylo_run_concord_extracted.log" --export=root=${root} ${root}/code/04_Phylogeny_and_Admixture/concord.sb
