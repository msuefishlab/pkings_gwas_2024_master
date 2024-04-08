## 02_Run_Prep_For_IQTree.sh
## Dec 2023
## JRG
## Submits to queue for Phylogeny Reconstruction

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

mkdir -p ${root}/output_data/slurm_logs/04_Phylogeny_And_Admixture/

alignments=("${root}/output_data/04_Phylogeny_And_Admixture/1000_rand_filtered_gene_alignments/"*.fasta)
nalignments=${#alignments[@]}

echo sbatch --job-name "PHYLO_IQTREE_LOCI" --output ${root}"/output_data/slurm_logs/04_Phylogeny_And_Admixture/phylo_run_iqtree_loci_%a.log" -a 0-$nalignments --export=root=${root} ${root}/code/04_Phylogeny_and_Admixture/run_iqtree_gene_loci_ind.sb
sbatch --job-name "PHYLO_IQTREE_LOCI" --output ${root}"/output_data/slurm_logs/04_Phylogeny_And_Admixture/phylo_run_iqtree_loci_%a.log" -a 0-$nalignments --export=root=${root} ${root}/code/04_Phylogeny_and_Admixture/run_iqtree_gene_loci_ind.sb
