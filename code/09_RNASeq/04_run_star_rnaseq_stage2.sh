mkdir -p ./stage2

for i in $SCRATCH/pkings_trimmed/*/tr_*1P*.fastq.gz
do
  sample=$(basename ${i//_1P.fastq.gz/});
  echo sbatch --output=${sample}-align-%j.out --export=READ1="${i}",READ2="${i//1P/2P}",OUTNAME=${sample} run_star_rnaseq_stage2.sb;
  sbatch --output=${sample}-align-%j.out --export=READ1="${i}",READ2="${i//1P/2P}",OUTNAME=${sample} run_star_rnaseq_stage2.sb;
 done
