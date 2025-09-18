#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source "$root/pkings_gwas.env"


slurmoutdir=${root}/output_data/slurm_outputs/09_RNASeq

mkdir -p $slurmoutdir

echo sbatch --output=${slurmoutdir}/STAR_INDEX_TRANSCRIPTOME.out ${root}/code/09_RNASeq/run_star_index.sb
sbatch --output=${slurmoutdir}/STAR_INDEX_TRANSCRIPTOME.out ${root}/code/09_RNASeq/run_star_index.sb
		