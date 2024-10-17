#!/bin/bash

## 01_create_region_lists.sh
## Feb 2023
## JRG
# This script creates region files for downstream phasing and is run interactively from the command line.

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

outdir=${root}/input_data/08_Phasing/
mkdir -p ${outdir}

mkdir -p ${outdir}/region_lists/


cat input_data/00_Reference_Genome/${reference}.fai | sort -r -k2 -n | head -25 | cut -f1 > ${outdir}/pseudo_chrs.txt

i=0
while read SCAFFOLD
do
  echo "$SCAFFOLD" > ${outdir}/region_lists/region_${i}.txt
  ((i++))
done < ${outdir}/pseudo_chrs.txt
cat input_data/00_Reference_Genome/${reference}.fai  | sort -r -k2 -n | tail -n +26 | cut -f1 > ${outdir}/region_lists/region_${i}.txt

