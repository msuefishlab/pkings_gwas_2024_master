#!/bin/bash

# 05_Prep_for_Astral.sh
# Aug 2024
# JRG
# Combines iqtree results for ASTRAL input

root="$(git rev-parse --show-toplevel)"
source "${root}/pkings_gwas.env"

cat ${root}/output_data/04_Phylogenomics/gene_trees_extracted/*treefile > ${root}/output_data/04_Phylogenomics/all_genes_extracted.treefile