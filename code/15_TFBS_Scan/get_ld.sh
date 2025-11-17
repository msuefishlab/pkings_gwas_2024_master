#!/bin/bash --login

root="$(git rev-parse --show-toplevel)"

source ${root}/pkings_gwas.env

outdir=/mnt/research/efish/pkings_gwas_2024_master/output_data/06_Association/PKINGS_ALL_WOB_EXCLUDED
OUTNAME=PKINGS_ALL_WOB_EXCLUDED

cat ${outdir}/${OUTNAME}_SNPS_IN_PEAKS.txt | awk '$13 >= 6 {print $0}' | cut -f2 | cut -d: -f1 > ${outdir}/${OUTNAME}_all_top_peak_snps.txt 


echo singularity exec --bind $root:/project_root --bind $outdir:/out_dir  ${gwas_tools_image} bash -c "bcftools annotate --threads 8 --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' -O z -o /project_root/sptbn4b_region_merge_ID.vcf.gz /project_root/sptbn4b_region.vcf.gz"
singularity exec --bind $root:/project_root --bind $outdir:/out_dir  ${gwas_tools_image} bash -c "bcftools annotate --threads 8 --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' -O z -o /project_root/sptbn4b_region_merge_ID.vcf.gz /project_root/sptbn4b_region.vcf.gz"

echo singularity exec --bind $root:/project_root --bind $outdir:/out_dir ${gwas_tools_image} bash -c "bcftools query -f '%ID\\:%CHROM\\:%POS, %POS\\, %CHROM\n' /project_root/sptbn4b_region_merge_ID.vcf.gz > /project_root/sptbn4b_region_merge_ID_merge_map.txt"
singularity exec --bind $root:/project_root --bind $outdir:/out_dir ${gwas_tools_image} bash -c "bcftools query -f '%ID\\:%CHROM\\:%POS, %POS\\, %CHROM\n' /project_root/sptbn4b_region_merge_ID.vcf.gz > /project_root/sptbn4b_region_merge_ID_merge_map.txt"


## First Make the PLINK Compatible Files
echo singularity exec --bind $root:/project_root --bind $outdir:/out_dir ${gwas_tools_image} /plink/plink \
--vcf /project_root/sptbn4b_region_merge_ID.vcf.gz \
--allow-extra-chr \
--allow-no-sex \
--chr-set 25 \
--out /project_root/sptbn4b_region_merge_ID.ld \
--r2 with-freqs \
--ld-window 10000 \
--ld-window-kb 10000 \
--ld-snp-list ${outdir}/${OUTNAME}_all_top_peak_snps.txt  \
--ld-window-r2 0.8 \
--threads 4

singularity exec --bind $root:/project_root --bind $outdir:/out_dir ${gwas_tools_image} /plink/plink \
--vcf /project_root/sptbn4b_region_merge_ID.vcf.gz \
--allow-extra-chr \
--allow-no-sex \
--chr-set 25 \
--out /project_root/sptbn4b_region_merge_ID.ld \
--r2 with-freqs \
--ld-window 10000 \
--ld-window-kb 10000 \
--ld-snp-list ${outdir}/${OUTNAME}_all_top_peak_snps.txt  \
--ld-window-r2 0.8 \
--threads 4
