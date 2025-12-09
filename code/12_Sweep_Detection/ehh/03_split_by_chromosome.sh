#!/bin/bash
set -euo pipefail

## split_vcf_by_chrom.sh
## Split polarized VCFs by chromosome for EHH analysis
##
## Usage:
##   bash code/12_Sweep_Detection/ehh/split_vcf_by_chrom.sh [--all-chromosomes]
##
## Options:
##   --all-chromosomes    Process all 25 chromosomes instead of just GWAS peak chromosomes
##
## Author: Jason Gallant Lab
## Date: 2024

# Parse command line arguments
ALL_CHROM=false
if [[ "${1:-}" == "--all-chromosomes" ]]; then
  ALL_CHROM=true
fi

module load BCFtools

root="$(git rev-parse --show-toplevel)"
outdir="${root}/output_data/12_Sweep_Detection"
vcfdir="${outdir}/vcfs"
outsub="${outdir}/by_chrom"

mkdir -p "${outsub}"
cd "${vcfdir}"

# Target chromosomes (with 'chr' prefix already)
if [[ "${ALL_CHROM}" == "true" ]]; then
  echo "Processing all 25 chromosomes"
  targets=(chr{1..25})
else
  echo "Processing GWAS peak chromosomes only (chr6, chr8, chr13, chr16, chr17, chr24)"
  targets=(chr6 chr8 chr13 chr16 chr17 chr24)
fi

shopt -s nullglob
vcfs=( *.polarized.vcf.gz )
shopt -u nullglob
[[ ${#vcfs[@]} -gt 0 ]] || { echo "No *.polarized.vcf.gz files found in ${vcfdir}"; exit 1; }

for vcf in "${vcfs[@]}"; do
  echo "Processing ${vcf} ..."
  base="${vcf%.vcf.gz}"   # e.g., BP1.polarized

  for region in "${targets[@]}"; do
    out="${outsub}/${base}.${region}.vcf.gz"
    echo "  - Subsetting ${region} → $(basename "${out}")"
    bcftools view -r "${region}" -Oz -o "${out}" "${vcf}"
    bcftools index -f "${out}"
  done
done

echo "Done. Per-chromosome VCFs are in: ${outsub}"
