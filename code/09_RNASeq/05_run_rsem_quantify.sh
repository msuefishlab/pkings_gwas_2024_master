mkdir -p ./RSEM

for i in stage1/*Aligned.toTranscriptome.out.bam
do
  sample=$(basename ${i//Aligned.toTranscriptome.out.bam/});
  echo sbatch --output=${sample}-RSEM-%j.out --export=OUTNAME=${sample} run_rsem.sb;
  sbatch --output=${sample}-RSEM-%j.out --export=OUTNAME=${sample} run_rsem.sb;
 done
