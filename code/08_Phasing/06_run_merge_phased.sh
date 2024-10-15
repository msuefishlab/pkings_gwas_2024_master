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

mkdir -p $WORK_D/slurm_outputs/merge_phased/
sbatch --job-name "MERGE_PHASED" --output $WORK_D/slurm_outputs/merge_phased/"MERGE_PHASED.slurm_%j_%a.log" --export=scriptdir=$scriptdir $scriptdir/merge_phased.sb