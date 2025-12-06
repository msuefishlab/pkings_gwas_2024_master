#!/bin/bash
set -euo pipefail

## 01_prepare_vcfs.sh
## Verify polarized VCFs exist for SweeD analysis
## SweeD reads VCF directly - no conversion needed!
##
## Usage:
##   bash code/12_Sweep_Detection/sweed/01_prepare_vcfs.sh [--all-chromosomes]
##
## Options:
##   --all-chromosomes    Check all 25 chromosomes instead of just GWAS peak chromosomes
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
outdir="${root}/output_data/12_Sweep_Detection/sweed"
mkdir -p "${outdir}/vcfs"

# Population sets to process
POPS=(BP1 TP1 BP2 TP2 BP3)

# Chromosomes to check
if [[ "${ALL_CHROM}" == "true" ]]; then
  echo "Checking all 25 chromosomes"
  CHROMOSOMES=(chr{1..25})
else
  echo "Checking GWAS peak chromosomes only (chr6, chr8, chr13, chr16, chr17, chr24)"
  CHROMOSOMES=(chr6 chr8 chr13 chr16 chr17 chr24)
fi

# Log file
log_file="${outdir}/vcf_validation.log"
echo "SweeD VCF Validation Report" > "${log_file}"
echo "Generated: $(date)" >> "${log_file}"
echo "========================================" >> "${log_file}"
echo "" >> "${log_file}"

# Manifest file for downstream processing
manifest="${outdir}/vcfs/vcf_manifest.txt"
echo -e "population\tchromosome\tvcf_path\thas_ancestral_allele" > "${manifest}"

# Counters
total_vcfs=0
found_vcfs=0
missing_vcfs=0
vcfs_with_aa=0
vcfs_without_aa=0

echo "Validating VCFs for SweeD analysis..."
echo ""

# Check each population × chromosome combination
for pop in "${POPS[@]}"; do
  echo "========================================" | tee -a "${log_file}"
  echo "Population: ${pop}" | tee -a "${log_file}"
  echo "========================================" | tee -a "${log_file}"

  for chr in "${CHROMOSOMES[@]}"; do
    total_vcfs=$((total_vcfs + 1))

    # Try per-chromosome VCF first
    vcf_by_chr="${root}/output_data/12_Sweep_Detection/vcfs/by_chrom/${pop}.polarized.${chr}.vcf.gz"

    # If per-chromosome doesn't exist, note that we'll use full VCF with filtering
    if [[ ! -f "${vcf_by_chr}" ]]; then
      # Check if full VCF exists
      vcf_full="${root}/output_data/12_Sweep_Detection/vcfs/${pop}.polarized.vcf.gz"

      if [[ ! -f "${vcf_full}" ]]; then
        echo "  [MISSING] ${chr}: No VCF found" | tee -a "${log_file}"
        echo -e "${pop}\t${chr}\tMISSING\tUNKNOWN" >> "${manifest}"
        missing_vcfs=$((missing_vcfs + 1))
        continue
      fi

      vcf="${vcf_full}"
      vcf_type="full (will filter to ${chr})"
    else
      vcf="${vcf_by_chr}"
      vcf_type="per-chromosome"
    fi

    found_vcfs=$((found_vcfs + 1))

    # Check if VCF has INFO/AA field (ancestral allele annotation)
    if bcftools view -h "${vcf}" 2>/dev/null | grep -q '##INFO=<ID=AA'; then
      has_aa="YES"
      vcfs_with_aa=$((vcfs_with_aa + 1))
      status="[OK]"
    else
      has_aa="NO"
      vcfs_without_aa=$((vcfs_without_aa + 1))
      status="[WARNING]"
    fi

    # Log result
    echo "  ${status} ${chr}: ${vcf_type}, AA field: ${has_aa}" | tee -a "${log_file}"

    # Add to manifest
    echo -e "${pop}\t${chr}\t${vcf}\t${has_aa}" >> "${manifest}"
  done

  echo "" | tee -a "${log_file}"
done

# Summary
echo "========================================" | tee -a "${log_file}"
echo "VALIDATION SUMMARY" | tee -a "${log_file}"
echo "========================================" | tee -a "${log_file}"
echo "Total VCFs expected: ${total_vcfs}" | tee -a "${log_file}"
echo "VCFs found: ${found_vcfs}" | tee -a "${log_file}"
echo "VCFs missing: ${missing_vcfs}" | tee -a "${log_file}"
echo "VCFs with AA field: ${vcfs_with_aa}" | tee -a "${log_file}"
echo "VCFs without AA field: ${vcfs_without_aa}" | tee -a "${log_file}"
echo "" | tee -a "${log_file}"

if [[ ${vcfs_without_aa} -gt 0 ]]; then
  echo "WARNING: ${vcfs_without_aa} VCFs are missing ancestral allele (INFO/AA) annotations!" | tee -a "${log_file}"
  echo "SweeD requires polarized VCFs. Run polarize_vcf.sh to add AA annotations." | tee -a "${log_file}"
  echo "" | tee -a "${log_file}"
fi

if [[ ${missing_vcfs} -gt 0 ]]; then
  echo "ERROR: ${missing_vcfs} VCFs are missing!" | tee -a "${log_file}"
  echo "Run split_vcfs.sh and split_vcf_by_chrom.sh to create population VCFs." | tee -a "${log_file}"
  echo "" | tee -a "${log_file}"
  exit 1
fi

echo "Manifest file created: ${manifest}" | tee -a "${log_file}"
echo "Log file: ${log_file}" | tee -a "${log_file}"
echo "" | tee -a "${log_file}"

if [[ ${vcfs_with_aa} -eq ${found_vcfs} ]]; then
  echo "SUCCESS: All VCFs are ready for SweeD analysis!" | tee -a "${log_file}"
  echo "" | tee -a "${log_file}"
  echo "Next step: Run 02_calculate_background_sfs.sh" | tee -a "${log_file}"
else
  echo "WARNING: Some VCFs lack AA annotations. SweeD will skip non-polarized sites." | tee -a "${log_file}"
  echo "" | tee -a "${log_file}"
  echo "Consider running polarize_vcf.sh before proceeding." | tee -a "${log_file}"
fi
