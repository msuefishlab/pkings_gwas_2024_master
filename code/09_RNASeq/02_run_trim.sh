
for i in /mnt/research/efish/lab_data/rnaseq/original_data/mormyriformes/pk_apa_bam_gabon2019/PK*R1*gz
do
  sample=$(basename ${i//_R1_001.fastq.gz/});
  mkdir -p $SCRATCH/pkings_trimmed/${sample}
  echo sbatch --output=${sample}-trim-%j.out --export=READ1="${i}",READ2="${i//R1/R2}",OUTNAME=$SCRATCH/pkings_trimmed/${sample},SAMPLE=${sample} run_trim.sb;
  sbatch --output=${sample}-trim-%j.out --export=READ1="${i}",READ2="${i//R1/R2}",OUTNAME=$SCRATCH/pkings_trimmed/${sample},SAMPLE=${sample} run_trim.sb;
 done
