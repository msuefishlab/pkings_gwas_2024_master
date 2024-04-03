## 03_Run_Prep_For_IQTree.sh
## Dec 2023
## JRG
## Submits to queue to Prepare Input Data for Phylogeny and ADMIXTURE Analysis

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

mkdir -p ${root}/output_data/slurm_logs/04_Phylogeny_And_Admixture/

outdir=${root}/output_data/04_Phylogeny_And_Admixture/
refdir=${root}/input_data/00_Reference_Genome/

mkdir -p ${outdir}
mkdir -p ${refdir}/split_gffs

awk -v outdir="${refdir}/split_gffs" '{print $0 >> outdir"/"$1".gff"}' ${refdir}/${gff_file}

gff_files=($(find ${refdir}/split_gffs -type f -name "*.gff" -exec realpath {} \;))
gff_file_count=${#gff_files[@]}

echo sbatch --job-name "SPLIT_VCFS" --output ${root}"/output_data/slurm_logs/04_Phylogeny_And_Admixture/phylo_split_gffs_%a.log" -a 0-${gff_file_count} --export=root=${root},gff_files=${gff_files} ${root}/code/04_Phylogeny_and_Admixture/split_vcf_to_genes.sb
sbatch --job-name "SPLIT_VCFS" --output ${root}"/output_data/slurm_logs/04_Phylogeny_And_Admixture/phylo_split_gffs_%a.log" -a 0-${gff_file_count} --export=root=${root},gff_files=${gff_files} ${root}/code/04_Phylogeny_and_Admixture/split_vcf_to_genes.sb
