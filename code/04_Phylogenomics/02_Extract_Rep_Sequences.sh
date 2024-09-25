## 03_Run_Extract_Rep_Sequences
## Sept 2024
## JRG
## Extracts Representitive Sequences for IQTree

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

indir=${root}/input_data/04_Phylogeny_And_Admixture/
outdir=${root}/output_data/04_Phylogeny_And_Admixture/

mkdir -p ${outdir}/extracted_seqs



# Loop over each file in the directory
for FILE in "${outdir}"/gene_alignments/*.fasta; do
	bn=$(basename "$FILE")
	if [ -e "${outdir}/extracted_seqs/${bn}" ]; then
	echo "File ${bn} exists, skipping!"
	else
    singularity exec --bind $root:/project_root --bind $outdir:/outdir --bind $indir:/indir  ${gwas_tools_image} /seqkit/seqkit grep -i -f /indir/representitives.txt /outdir/gene_alignments/${bn} > ${outdir}/extracted_seqs/${bn}
	fi
done

