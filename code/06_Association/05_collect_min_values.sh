#!/bin/bash

## 04_collect_min_values.sh
## Feb 2023
## JRG
# This script collects minimum p-values across permuted gemma gwas runs, which would normally be collected by gemma-wrapper when not in slurm mode.
# e.g. . code/06_Association/05_collect_min_values.sh APA_BAM_ONLY_BP1_TP0_WOBBLE9

root="$(git rev-parse --show-toplevel)"


OUTNAME=$1

cd $root

source $root/pkings_gwas.env

indir=${root}/output_data/06_Association/
outdir=${root}/output_data/06_Association/$OUTNAME/${OUTNAME}_permution/

rm -rf ${outdir}/min.txt
touch ${outdir}/min.txt
for f in ${outdir}/*.assoc.txt; 
	do  
		if grep -e "phenotypes-" ${f%%.assoc.txt}.log.txt; 
			then echo "Processing $f file..."; 
				awk 'NR==1 || $10 < min { line = $0; min = $10}; END {print line}' $f >> ${outdir}/min.txt; 
		fi; 
	done