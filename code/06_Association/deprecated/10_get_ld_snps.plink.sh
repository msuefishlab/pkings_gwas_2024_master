#!/bin/bash --login

OUTNAME=$1

root="$(git rev-parse --show-toplevel)"

source $root/pkings_gwas.env

outdir=${root}/output_data/06_Association/${OUTNAME}

mkdir -p ${outdir}

cat ${outdir}/${OUTNAME}_PLINK_ALLELIC_SNPS_IN_PEAKS.txt | cut -f2 | cut -d: -f1 > ${outdir}/${OUTNAME}_all_plink_allelic_peak_snps.txt 

## First Make the PLINK Compatible Files
echo singularity exec --bind $root:/project_root --bind $outdir:/out_dir ${gwas_tools_image} /plink/plink \
--vcf /out_dir/${OUTNAME}_merge_ID.vcf.gz \
--allow-extra-chr \
--allow-no-sex \
--chr-set 25 \
--out /out_dir/${OUTNAME}_plink_allelic_index_snps.ld \
--r2 with-freqs \
--ld-window 10000 \
--ld-window-kb 10000 \
--ld-snp-list ${outdir}/${OUTNAME}_all_plink_allelic_peak_snps.txt \
--ld-window-r2 0.25 \
--threads 4

singularity exec --bind $root:/project_root --bind $outdir:/out_dir ${gwas_tools_image} /plink/plink \
--vcf /out_dir/${OUTNAME}_merge_ID.vcf.gz \
--allow-extra-chr \
--allow-no-sex \
--chr-set 25 \
--out /out_dir/${OUTNAME}_plink_allelic_index_snps.ld \
--r2 with-freqs \
--ld-window 10000 \
--ld-window-kb 10000 \
--ld-snp-list ${outdir}/${OUTNAME}_all_plink_allelic_peak_snps.txt \
--ld-window-r2 0.25 \
--threads 4