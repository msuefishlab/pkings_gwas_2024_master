## run_merge.sh
## 12-03-23
## JRG
## Merges interval split tab delimited files for QC Analysis

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

outdir=${scratch_store}/output_data/03_QC/

mkdir -p ${outdir}

awk 'FNR==1 && NR!=1{next;}{print}' ${outdir}/*indel*tab > ${root}/output_data/03_QC/${vcf_file%.vcf.gz}.all.indel.tab
awk 'FNR==1 && NR!=1{next;}{print}' ${outdir}/*snp*tab > ${root}/output_data/03_QC/${vcf_file%.vcf.gz}.all.snp.tab