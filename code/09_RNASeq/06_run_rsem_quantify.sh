#!/bin/bash

## Run by bash code/09_RNASeq/02_cat_by_specno.sh input_data/09_RNASeq/rnaseq_metadata.txt

root="$(git rev-parse --show-toplevel)"
source "$root/pkings_gwas.env"

# Define the input file
input_file=$1

indir="${root}/output_data/09_RNASeq/aligned_reads"  # Define the output director
outdir="${root}/output_data/09_RNASeq/quantification"  # Define the output director
slurmoutdir=${root}/output_data/slurm_outputs/09_RNASeq

mkdir -p $outdir
mkdir -p $slurmoutdir

# Skip the header and read the file line by line using a file descriptor
{
  read -r header # Skip the header
  declare -A seen_specnos # Declare an associative array to track unique specnos

  while IFS=$'\t' read -r specno population r1 r2 phenotype; do
    if [[ -z "${seen_specnos[$specno]}" ]]; then
      # If specno is not yet seen
      seen_specnos[$specno]=1 # Mark specno as seen
    echo sbatch --output=${slurmoutdir}/RSEM_${specno}.out --export=INDIR=${indir}/${specno},OUTDIR=${outdir}/${specno},SAMPLE=${specno} ${root}/code/09_RNASeq/run_rsem.sb;
    sbatch --output=${slurmoutdir}/RSEM_${specno}.out --export=INDIR=${indir}/${specno},OUTDIR=${outdir}/${specno},SAMPLE=${specno} ${root}/code/09_RNASeq/run_rsem.sb;
   fi
  done
} < "$input_file"
