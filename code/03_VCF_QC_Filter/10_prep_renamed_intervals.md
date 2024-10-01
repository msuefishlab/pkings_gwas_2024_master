outdir=${root}/input_data/03_QC
indir=$root/input_data/00_Reference_Genome/

mkdir -p ${outdir}

source pkings_gwas.env
docker run -d -it --name gatk --mount type=bind,source=$root/input_data/00_Reference_Genome/,target=/input_data --mount type=bind,source=$root/input_data/03_QC/,target=/output_data broadinstitute/gatk


singularity shell --bind $root:/project_root --bind $outdir:/output_data --bind $indir:/input_data ${gatk_image}

#### prepare files/

cat /input_data/$reference.fai | cut -f1 > /output_data/intervals.renamed.list

/gatk/gatk CreateSequenceDictionary R= /input_data/$reference

/gatk/gatk SplitIntervals -R /input_data/$reference -L /output_data/intervals.renamed.list -scatter 50 -O /output_data/interval-files -mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW

cat /input_data/$reference.fai | cut -f1,2 | sort -nrk 2 | head -25 | cut -f1 > /output_data/chrom_intervals.list
cat /input_data/$reference.fai | cut -f1,2 | sort -nrk 2 | tail +26 | cut -f1 > /output_data/remaining_intervals.list

/gatk/gatk SplitIntervals -R /input_data/$reference -L /output_data/chrom_intervals.list --scatter-count 100 -O /output_data/chrom-interval-files --subdivision-mode INTERVAL_SUBDIVISION

/gatk/gatk SplitIntervals -R /input_data/$reference -L /output_data/remaining_intervals.list --scatter-count 100 -O /output_data/remaining-interval-files --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW