## 03_Run_Prep_For_IQTree.sh
## Dec 2023
## JRG
## Submits to queue to Prepare Input Data for Phylogeny and ADMIXTURE Analysis

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

mkdir -p ${root}/output_data/slurm_logs/04_Phylogeny_And_Admixture/

echo sbatch --job-name "PREP_FOR_PHYLO_ADMIX" --output ${root}"/output_data/slurm_logs/04_Phylogeny_And_Admixture/phylo_prep_iq_tree.log" --export=root=${root} ${root}/code/04_Phylogeny_And_Admixture/prep_for_iqtree_and_admixture.sb
sbatch --job-name "PREP_FOR_PHYLO_ADMIX" --output ${root}"/output_data/slurm_logs/04_Phylogeny_And_Admixture/phylo_prep_iq_tree.log" --export=root=${root} ${root}/code/04_Phylogeny_And_Admixture/prep_for_iqtree_and_admixture.sb
