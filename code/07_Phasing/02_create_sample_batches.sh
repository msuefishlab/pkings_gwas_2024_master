#!/bin/bash

## 01_create_region_lists.sh
## Feb 2023
## JRG
# This script creates region files for downstream phasing and is run interactively from the command line.

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

outdir=${root}/input_data/07_Phasing/
mkdir -p ${outdir}

module purge; module load BCFtools

bcftools query -l output_data/03_QC/${snps_only_for_phasing} | split -l 25 -d - input_data/07_Phasing/sample_names_b

