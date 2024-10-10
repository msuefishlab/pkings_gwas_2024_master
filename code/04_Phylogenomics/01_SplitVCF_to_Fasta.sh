## 03_Run_Prep_For_IQTree.sh
## Dec 2023
## JRG
## Submits to queue to Prepare Input Data for Phylogeny and ADMIXTURE Analysis

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

mkdir -p ${root}/output_data/slurm_logs/04_Phylogenomics/

outdir=${root}/output_data/04_Phylogenomics/
refdir=${root}/input_data/00_Reference_Genome/

mkdir -p ${outdir}
mkdir -p ${refdir}/split_regions

# Calculate lines per file for approximately 200 files
total_lines=$(wc -l < ${refdir}/pkings_gene_regions.txt)
lines_per_file=$((total_lines / 200 + 1))

# Split the file with numeric suffixes
split -d -l "$lines_per_file" ${refdir}/pkings_gene_regions.txt ${refdir}/split_regions/chunk_ -a 3 --numeric-suffixes=1

# This will create files named chunk_001, chunk_002, ..., each with around the same number of lines

echo sbatch --job-name "SPLIT_VCFS" --output ${root}"/output_data/slurm_logs/04_Phylogenomics/phylo_split_gffs_%a.log" -a 1-200 --export=root=${root} ${root}/code/04_Phylogenomics//split_vcf_to_genes.sb
sbatch --job-name "SPLIT_VCFS" --output ${root}"/output_data/slurm_logs/04_Phylogenomics/phylo_split_gffs_%a.log" -a 1-200 --export=root=${root} ${root}/code/04_Phylogenomics//split_vcf_to_genes.sb
