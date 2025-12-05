#!/bin/bash
set -euo pipefail

## 01_convert_vcf_to_sf2.sh
## Convert polarized VCFs to SweepFinder2 format
## Per population, per chromosome
##
## Usage:
##   bash code/12_Sweep_Detection/01_convert_vcf_to_sf2.sh [--all-chromosomes]
##
## Options:
##   --all-chromosomes    Process all 25 chromosomes instead of just GWAS peak chromosomes
##
## Author: Jason Gallant Lab
## Date: 2024

root="$(git rev-parse --show-toplevel)"
source "${root}/pkings_gwas.env"

# Parse command line arguments
ALL_CHROM=false
if [[ "${1:-}" == "--all-chromosomes" ]]; then
  ALL_CHROM=true
fi

# Create output directories
outdir="${root}/output_data/12_Sweep_Detection/sweepfinder2"
mkdir -p "${outdir}/sfs_input"

# Population sets to process (actual populations, not comparison groups)
POPS=(BP1 TP1 BP2 TP2 BP3)

# Chromosomes to process
if [[ "${ALL_CHROM}" == "true" ]]; then
  echo "Processing all 25 chromosomes"
  CHROMOSOMES=(chr{1..25})
else
  echo "Processing GWAS peak chromosomes only (chr6, chr8, chr13, chr16, chr17, chr24)"
  CHROMOSOMES=(chr6 chr8 chr13 chr16 chr17 chr24)
fi

# Log file
log_file="${outdir}/conversion.log"
echo "Starting VCF to SweepFinder2 conversion: $(date)" > "${log_file}"
echo "Processing ${#POPS[@]} populations × ${#CHROMOSOMES[@]} chromosomes = $((${#POPS[@]} * ${#CHROMOSOMES[@]})) files" >> "${log_file}"

# Counter for progress
total_files=$((${#POPS[@]} * ${#CHROMOSOMES[@]}))
current=0

# Process each population
for pop in "${POPS[@]}"; do
  echo ""
  echo "========================================" | tee -a "${log_file}"
  echo "Processing population: ${pop}" | tee -a "${log_file}"
  echo "========================================" | tee -a "${log_file}"

  for chr in "${CHROMOSOMES[@]}"; do
    current=$((current + 1))
    echo ""
    echo "[$current/$total_files] Converting ${pop} ${chr}..." | tee -a "${log_file}"

    # Input VCF path (per-chromosome subset)
    vcf_by_chr="${root}/output_data/12_Sweep_Detection/by_chrom/${pop}.polarized.${chr}.vcf.gz"

    # If per-chromosome VCF doesn't exist, try to create it from full VCF
    if [[ ! -f "${vcf_by_chr}" ]]; then
      echo "  WARNING: ${vcf_by_chr} not found" | tee -a "${log_file}"

      # Try full VCF
      vcf_full="${root}/output_data/12_Sweep_Detection/${pop}.polarized.vcf.gz"

      if [[ ! -f "${vcf_full}" ]]; then
        echo "  ERROR: Neither per-chromosome nor full VCF found for ${pop}" | tee -a "${log_file}"
        echo "    Tried: ${vcf_by_chr}" | tee -a "${log_file}"
        echo "    Tried: ${vcf_full}" | tee -a "${log_file}"
        continue
      fi

      echo "  Using full VCF and filtering to ${chr}" | tee -a "${log_file}"
      vcf="${vcf_full}"
      chr_arg="--chromosome ${chr}"
    else
      vcf="${vcf_by_chr}"
      chr_arg=""
    fi

    # Output SweepFinder2 format file
    output="${outdir}/sfs_input/${pop}.${chr}.sf2.txt"

    # Run Python conversion script
    echo "  Input:  ${vcf}" | tee -a "${log_file}"
    echo "  Output: ${output}" | tee -a "${log_file}"

    python3 "${root}/code/12_Sweep_Detection/vcf_to_sweepfinder2.py" \
      "${vcf}" "${output}" ${chr_arg} 2>&1 | tee -a "${log_file}"

    # Check if output was created and has content
    if [[ -f "${output}" ]]; then
      line_count=$(wc -l < "${output}")
      echo "  Output file has ${line_count} lines (including header)" | tee -a "${log_file}"
    else
      echo "  ERROR: Output file not created!" | tee -a "${log_file}"
    fi
  done
done

echo "" | tee -a "${log_file}"
echo "========================================" | tee -a "${log_file}"
echo "Conversion complete: $(date)" | tee -a "${log_file}"
echo "========================================" | tee -a "${log_file}"
echo "" | tee -a "${log_file}"
echo "Output directory: ${outdir}/sfs_input/" | tee -a "${log_file}"
echo "Log file: ${log_file}" | tee -a "${log_file}"
echo "" | tee -a "${log_file}"
echo "Files created:" | tee -a "${log_file}"
ls -lh "${outdir}/sfs_input/"*.sf2.txt | wc -l | tee -a "${log_file}"
echo "" | tee -a "${log_file}"

# Summary statistics
echo "Summary by population:" | tee -a "${log_file}"
for pop in "${POPS[@]}"; do
  count=$(ls "${outdir}/sfs_input/${pop}.chr"*.sf2.txt 2>/dev/null | wc -l)
  echo "  ${pop}: ${count} chromosomes" | tee -a "${log_file}"
done

echo "" | tee -a "${log_file}"
echo "Next step: Run 02_calculate_background_sfs.sh to calculate genome-wide background SFS" | tee -a "${log_file}"
