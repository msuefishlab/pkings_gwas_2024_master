#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

mkdir -p ${root}/output_data/slurm_logs/09_RNASeq/

outdir=${root}/output_data/09_RNASeq/
indir=${root}/input_data/09_RNASeq/

mkdir -p ${outdir}

echo -n "" > ${outdir}/missing_files.txt
sed 1d ${indir}/sample_table_for_kallisto.csv | while read line; do

  export individual=$(echo ${line} | cut -d, -f1)
  export species=$(echo ${line} | cut -d, -f2)
  export tissue=$(echo ${line} | cut -d, -f3)
  export transcriptome=$(echo ${line} | cut -d, -f4)
  export reads=$(echo ${line} | cut -d, -f5)

  reads_array=($reads)

  if [ -f ${reads_array[1]} ]; then
      echo sbatch --job-name=$species.$tissue.$individual --output "${root}/output_data/slurm_logs/09_RNASeq/$species.$tissue.$individual.log" --export=individual,species,tissue,reads,transcriptome,root ${root}/code/09_RNASeq/submit_kallisto.sb
      sbatch --job-name=$species.$tissue.$individual --output "${root}/output_data/slurm_logs/09_RNASeq/$species.$tissue.$individual.log" --export=individual,species,tissue,reads,transcriptome,root ${root}/code/09_RNASeq/submit_kallisto.sb
    else
      echo "${reads_array[1]} does not exist."
      echo ${reads} >> missing_files.txt
  fi

done
