#!/bin/bash

## get_SNP_HITS.sh
## February 2024
## JRG
# This script creates a genotype matrix for a top100 snps file generated by GWAS report pipeline:
# bash 01_get_SNP_HITS.sh $SCRATCH/pkings_popgen_2023/output_data/gemma/apa_bam_no_wobbles/peak_data/apa_bam_no_wobbles_OUTLIER_SNPs_TOP100.txt


root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

OUTNAME=$1

outdir=${root}/output_data/08_Peak_Analysis/${OUTNAME}
mkdir -p ${outdir}



TOP_SNPS=${root}/output_data/06_Association/${OUTNAME}/${OUTNAME}"_SNPS_IN_PEAKS_TOP250.txt"


cut -f 1,2 $TOP_SNPS > ${outdir}/${OUTNAME}_snps.txt

singularity exec ${gwas_tools_image} bcftools view -R ${outdir}/${OUTNAME}_snps.txt ${root}/output_data/03_QC/$snps_only_vcf -O z -o ${outdir}/${OUTNAME}_SNPS_IN_PEAKS_TOP250.vcf.gz

singularity exec ${gwas_tools_image} bash -c "bcftools query -l ${outdir}/${OUTNAME}_SNPS_IN_PEAKS_TOP250.vcf.gz | tr \"\\n\" \"\\t\" > ${outdir}/${OUTNAME}_header.txt"

echo "" >>  ${outdir}/${OUTNAME}_header.txt #add new line at end of header

singularity exec ${gwas_tools_image} bash -c "bcftools query -f '[%GT\t]\n' ${outdir}/${OUTNAME}_SNPS_IN_PEAKS_TOP250.vcf.gz > ${outdir}/${OUTNAME}_raw.genos"
cat ${outdir}/${OUTNAME}_header.txt ${outdir}/${OUTNAME}_raw.genos >${outdir}/${OUTNAME}.genos

singularity exec ${gwas_tools_image} bash -c "bcftools query -f '%CHROM\t %POS\t %REF\t %ALT\n' ${outdir}/${OUTNAME}_SNPS_IN_PEAKS_TOP250.vcf.gz > ${outdir}/${OUTNAME}_raw.pos"

echo -e "CHROM\tPOS\tREF\tALT" | cat - ${outdir}/${OUTNAME}_raw.pos > ${outdir}/${OUTNAME}.pos

paste -d "\t" ${outdir}/${OUTNAME}.pos ${outdir}/${OUTNAME}.genos > ${outdir}/${OUTNAME}.pos.genos

sed ${outdir}/${OUTNAME}.pos.genos \
  -e "s:0|1:1:g" \
  -e "s:1|0:1:g" \
  -e "s:0|0:0:g" \
  -e "s:1|1:2:g" \
  -e "s:0/1:1:g" \
  -e "s:1/0:1:g" \
  -e "s:0/0:0:g" \
  -e "s:1/1:2:g" \
  -e "s:./.:NA:g" \
  -e "s:.|.:NA:g" > ${outdir}/${OUTNAME}.bin.pos.genos

rm ${outdir}/${OUTNAME}.genos
rm ${outdir}/${OUTNAME}.pos
rm ${outdir}/${OUTNAME}_header.txt
rm ${outdir}/${OUTNAME}*_raw*
