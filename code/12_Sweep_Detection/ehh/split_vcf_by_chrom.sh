#!/bin/bash
set -euo pipefail

## split_each_pop_by_chrom.sh (simple)
## JRG — August 2025

module load BCFtools

root="$(git rev-parse --show-toplevel)"
outdir="${root}/output_data/12_Sweep_Detection"
vcfdir="${outdir}/vcfs"
outsub="${vcfdir}/by_chrom"

mkdir -p "${outsub}"
cd "${vcfdir}"

# target chromosomes (with 'chr' prefix already)
targets=(chr24 chr17 chr6 chr16 chr13 chr8)

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
