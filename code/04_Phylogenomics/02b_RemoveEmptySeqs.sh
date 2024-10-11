## 03_Run_Extract_Rep_Sequences
## Sept 2024
## JRG
## Extracts Representitive Sequences for IQTree

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

indir=${root}/input_data/04_Phylogenomics/
outdir=${root}/output_data/04_Phylogenomics/

mkdir -p ${outdir}/extracted_complete_alignments



find ${root}/output_data/04_Phylogenomics/extracted_seqs/ -type f -size +0c -exec cp {} ${root}/output_data/04_Phylogenomics/extracted_complete_alignments/ \;

