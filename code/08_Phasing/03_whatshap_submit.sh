#!/bin/bash

## 02_whatshap_submit.sh
## Feb 2023
## JRG
# This script creates runs whatshap on sample batches located in the metadata folder
# e.g. bash 03_whatshap_submit.sh sample_names_b_0.txt

scriptdir="$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )"
WORK_D=$(git rev-parse --show-toplevel)

cd $WORK_D

echo $scriptdir

source $WORK_D/code/00_utility/parse_yaml.sh

eval $(parse_yaml ./global_params.yaml)

while read SAMPLE
do
  cd $(git rev-parse --show-toplevel)
  echo "$SAMPLE"
  mkdir -p ./output_data/whatshap/$SAMPLE
  cd ./output_data/whatshap/$SAMPLE
  mkdir -p $WORK_D/slurm_outputs/whatshap
  sbatch --job-name $SAMPLE"_WHATSHAP" --array=0-25 --output $WORK_D/slurm_outputs/whatshap/$SAMPLE"_WHATSHAP.slurm_%j_%a.log" --export=scriptdir=$scriptdir $scriptdir/submit_whatshap.sb $SAMPLE

done < "$1"