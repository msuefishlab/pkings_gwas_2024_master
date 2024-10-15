#!/bin/bash

## 02_whatshap_submit.sh
## Feb 2023
## JRG
# This script runs shapeit on all whatshap vcf files
# e.g. bash 05_run_shapeit.sh

scriptdir="$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )"
WORK_D=$(git rev-parse --show-toplevel)

cd $WORK_D

echo $scriptdir

source $WORK_D/code/00_utility/parse_yaml.sh

eval $(parse_yaml ./global_params.yaml)

cp ./metadata/pseudo_chrs.txt ./metadata/shapeit_regions.txt
echo "final_chunk" >> ./metadata/shapeit_regions.txt

num_regions=$(cat ./metadata/shapeit_regions.txt | wc -l)


mkdir -p $WORK_D/slurm_outputs/shapeit/
sbatch --job-name $SAMPLE"_SHAPEIT" -a 1-${num_regions} --output $WORK_D/slurm_outputs/shapeit/"SHAPEIT.slurm_%j_%a.log" --export=scriptdir=$scriptdir $scriptdir/run_shapeit.sb