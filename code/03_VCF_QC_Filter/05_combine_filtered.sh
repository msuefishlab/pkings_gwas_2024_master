#!/bin/bash

## 05_combine_filtered.sh
## Dec 2023
## JRG
# This script submits bcftools concat jobs to slurm  to create downstream VCF files for analysis.

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"
scratchdir=${scratch_store}/output_data/03_QC/
outdir=${root}/output_data/03_QC/

ls -1 ${scratchdir}*.SNP.ONLY.BIALLELIC.ONLY.filtered.setGT.PASS.vcf.gz > ${outdir}/biallelic_only_shards.txt
ls -1 ${scratchdir}*.SNP.BIALLELIC.ONLY.AND.INV.TOGETHER.filtered.setGT.PASS.vcf.gz > ${outdir}/biallelic_and_invariant_shards.txt
ls -1 ${scratchdir}*.INDEL.filtered.setGT.vcf.gz > ${outdir}/indel_shards.txt

cat ${outdir}/biallelic_only_shards.txt > ${outdir}/FILTERED_BIALLELIC_SNPS_ONLY.list
cat ${outdir}/biallelic_only_shards.txt ${outdir}/indel_shards.txt > ${outdir}/FILTERED_BIALLELIC_AND_INDEL.list
cat ${outdir}/biallelic_and_invariant_shards.txt > ${outdir}/FILTERED_BIALLELIC_AND_INVARIANT.list

echo sbatch --job-name "COMBINE_FILTERED" --output ${root}"/output_data/slurm_logs"/03_QC/"qc_combined_filtered_%a.log " -a 0-2 --export=root=${root} ${root}/code/03_VCF_QC_Filter/combine_filtered.sb
sbatch --job-name "COMBINE_FILTERED" --output ${root}"/output_data/slurm_logs"/03_QC/"qc_combined_filtered_%a.log " -a 0-2 --export=root=${root} ${root}/code/03_VCF_QC_Filter/combine_filtered.sb
