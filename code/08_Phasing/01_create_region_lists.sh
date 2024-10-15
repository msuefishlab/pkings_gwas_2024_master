#!/bin/bash

## 01_create_region_lists.sh
## Feb 2023
## JRG
# This script creates region files for downstream phasing and is run interactively from the command line.

scriptdir="$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )"
WORK_D=$(git rev-parse --show-toplevel)

cd $WORK_D

echo $scriptdir

source $WORK_D/code/00_utility/parse_yaml.sh

eval $(parse_yaml ./global_params.yaml)


mkdir -p ./metadata/region_lists

cat ./data/reference/${REFERENCE_FILE}.fai | sort -r -k2 -n | head -25 | cut -f1 > ./metadata/pseudo_chrs.txt

i=0
while read SCAFFOLD
do
  echo "$SCAFFOLD" > ./metadata/region_lists/region_${i}.txt
  ((i++))
done < ./metadata/pseudo_chrs.txt
cat data/reference/ragtag.scaffold.fasta.fai | sort -r -k2 -n | tail -n +26 | cut -f1 > ./metadata/region_lists/region_${i}.txt

