#!/bin/bash

## 02_whatshap_submit.sh
## Feb 2023
## JRG
# This script combines all unplaced scaffolds into a single VCF to reduce the total number of files, and improve slurm throughput for shapeit
# e.g. bash 04_squish_outputs.sh

scriptdir="$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )"
WORK_D=$(git rev-parse --show-toplevel)

cd $WORK_D

echo $scriptdir

source $WORK_D/code/00_utility/parse_yaml.sh

eval $(parse_yaml ./global_params.yaml)

#cat $WORK_D/metadata/sample_names_b*.txt > $WORK_D/metadata/all_samples.txt

while read SAMPLE
do
  echo "$SAMPLE"
  mkdir -p $WORK_D/slurm_outputs/squish/
  sbatch --job-name $SAMPLE"_SQUISH" --output $WORK_D/slurm_outputs/squish/$SAMPLE"_SQUISH.slurm_%j_%a.log" --export=scriptdir=$scriptdir $scriptdir/squish_outputs.sb $SAMPLE

done < $WORK_D/metadata/stragglers.txt