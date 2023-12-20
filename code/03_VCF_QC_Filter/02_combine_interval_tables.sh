## run_merge.sh
## 12-03-23
## JRG
## Merges interval split tab delimited files for QC Analysis

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

awk 'FNR==1 && NR!=1{next;}{print}' ${root}/output_data/03_QC/*indel*tab > ${root}/output_data/03_QC/${vcf_file%.vcf.gz}.all.indel.tab
awk 'FNR==1 && NR!=1{next;}{print}' ${root}/output_data/03_QC/*snp*tab > ${root}/output_data/03_QC/${vcf_file%.vcf.gz}.all.snp.tab