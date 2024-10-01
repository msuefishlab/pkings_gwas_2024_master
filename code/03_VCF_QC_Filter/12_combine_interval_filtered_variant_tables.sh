## 02_combine_internal_tables.sh
## 12-03-23
## JRG
## Combines variant tables into one large table for QC analysis

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

outdir=${scratch_store}/output_data/03_QC/

mkdir -p ${outdir}

awk 'FNR==1 && NR!=1{next;}{print}' ${outdir}/${vcf_file%.vcf.gz}.*.for_inspection.tab > ${root}/output_data/03_QC/${snps_only_vcf%.vcf.gz}.filtered.stats.tab