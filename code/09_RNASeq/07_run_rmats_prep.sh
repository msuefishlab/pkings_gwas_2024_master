#!/bin/bash

## Run by bash code/09_RNASeq/02_cat_by_specno.sh input_data/09_RNASeq/rnaseq_metadata.txt

root="$(git rev-parse --show-toplevel)"
source "$root/pkings_gwas.env"

# Define the input file
input_file=$1

bamdir="${root}/output_data/09_RNASeq/aligned_reads_genome"  # Define the output director
outdir="${root}/output_data/09_RNASeq/rmats"  # Define the output director
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
    mkdir -p ${outdir}/${specno};
	bamfile=${bamdir}/${specno}/${specno}Aligned.sortedByCoord.out.bam
	echo $bamfile > ${outdir}/${specno}/bamlist.txt
    echo sbatch --output=${slurmoutdir}/RMATS_prep_${specno}.out --export=BAMPATH="${outdir}/${specno}/bamlist.txt",OUTDIR=${outdir},SAMPLE=${specno} ${root}/code/09_RNASeq/run_rmats_prep.sb;
	sbatch --output=${slurmoutdir}/RMATS_prep_${specno}.out --export=BAMPATH="${outdir}/${specno}/bamlist.txt",OUTDIR=${outdir},SAMPLE=${specno} ${root}/code/09_RNASeq/run_rmats_prep.sb;
 fi
  done
} < "$input_file"


