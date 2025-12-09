#!/bin/bash
set -euo pipefail

## 03_merge_ehh_results.sh
## Merge per-chromosome inES/iHS and Rsb results into genome-wide files
##
## Usage:
##   bash code/12_Sweep_Detection/ehh/03_merge_ehh_results.sh
##
## Prerequisites:
##   - 01_submit_ines_scans.sh must be completed (inES files must exist)
##   - 02_submit_rsb_scans.sh must be completed (Rsb files must exist)
##
## Author: Jason Gallant Lab
## Date: 2024

root="$(git rev-parse --show-toplevel)"
source "${root}/pkings_gwas.env"

# Define directories
ines_dir="${root}/output_data/12_Sweep_Detection/ehh/ines_scans"
rsb_dir="${root}/output_data/12_Sweep_Detection/ehh/rsb_scans"
outdir="${root}/output_data/12_Sweep_Detection/ehh/merged"
mkdir -p "${outdir}"

# Population sets
POPS=(BP1 BP2 BP3 TP1 TP2)

# Population comparisons for Rsb
COMPARISONS=(
  "BP1_TP1"
  "BP2_TP2"
  "BP3_TP2"
)

# Log file
log_file="${outdir}/merge.log"
echo "EHH Results Merging" > "${log_file}"
echo "Generated: $(date)" >> "${log_file}"
echo "========================================" >> "${log_file}"
echo "" >> "${log_file}"

echo "========================================"
echo "Merging per-chromosome EHH results"
echo "========================================"
echo ""

# ========================================
# Part 1: Merge inES/iHS files
# ========================================

echo "Part 1: Merging inES/iHS results..." | tee -a "${log_file}"
echo "" | tee -a "${log_file}"

for pop in "${POPS[@]}"; do
  echo "========================================" | tee -a "${log_file}"
  echo "Population: ${pop}" | tee -a "${log_file}"
  echo "========================================" | tee -a "${log_file}"

  # Output merged file
  merged_file="${outdir}/${pop}.genome_wide_ines.txt.gz"

  # Find all chromosome files for this population
  chr_files=$(ls ${ines_dir}/${pop}.chr*.ines.txt.gz 2>/dev/null || true)

  if [[ -z "${chr_files}" ]]; then
    echo "  WARNING: No inES files found for ${pop}" | tee -a "${log_file}"
    continue
  fi

  # Count files
  file_count=$(echo "${chr_files}" | wc -w)
  echo "  Found ${file_count} chromosome files" | tee -a "${log_file}"

  # Create temporary file for merging
  temp_file="${outdir}/.temp_${pop}_ines.txt"

  # Process first file to get header
  first_file=$(echo "${chr_files}" | awk '{print $1}')
  zcat "${first_file}" | head -n 1 > "${temp_file}"

  # Merge all chromosomes (skip headers for all files)
  total_lines=0

  for chr_file in ${chr_files}; do
    chr=$(basename "${chr_file}" .ines.txt.gz | sed "s/${pop}\.//")
    echo "    Processing ${chr}..." | tee -a "${log_file}"

    # Extract data lines (skip header)
    line_count=$(zcat "${chr_file}" | tail -n +2 | tee -a "${temp_file}" | wc -l)
    total_lines=$((total_lines + line_count))

    echo "      Added ${line_count} lines" | tee -a "${log_file}"
  done

  # Compress merged file
  gzip -c "${temp_file}" > "${merged_file}"
  rm -f "${temp_file}"

  # Summary for this population
  echo "" | tee -a "${log_file}"
  echo "  Created: ${merged_file}" | tee -a "${log_file}"
  echo "  Total data lines: ${total_lines}" | tee -a "${log_file}"
  echo "  File size: $(du -h ${merged_file} | cut -f1)" | tee -a "${log_file}"

  # Show first few lines
  echo "" | tee -a "${log_file}"
  echo "  First 5 data lines:" | tee -a "${log_file}"
  zcat "${merged_file}" | head -n 6 | tail -n 5 | tee -a "${log_file}"

  echo "" | tee -a "${log_file}"
done

# ========================================
# Part 2: Merge Rsb files
# ========================================

echo "" | tee -a "${log_file}"
echo "Part 2: Merging Rsb results..." | tee -a "${log_file}"
echo "" | tee -a "${log_file}"

for comp in "${COMPARISONS[@]}"; do
  echo "========================================" | tee -a "${log_file}"
  echo "Comparison: ${comp}" | tee -a "${log_file}"
  echo "========================================" | tee -a "${log_file}"

  # Output merged file
  merged_file="${outdir}/${comp}.genome_wide_rsb.txt.gz"

  # Find all chromosome files for this comparison
  chr_files=$(ls ${rsb_dir}/${comp}.chr*.rsb.txt.gz 2>/dev/null || true)

  if [[ -z "${chr_files}" ]]; then
    echo "  WARNING: No Rsb files found for ${comp}" | tee -a "${log_file}"
    continue
  fi

  # Count files
  file_count=$(echo "${chr_files}" | wc -w)
  echo "  Found ${file_count} chromosome files" | tee -a "${log_file}"

  # Create temporary file for merging
  temp_file="${outdir}/.temp_${comp}_rsb.txt"

  # Process first file to get header
  first_file=$(echo "${chr_files}" | awk '{print $1}')
  zcat "${first_file}" | head -n 1 > "${temp_file}"

  # Merge all chromosomes (skip headers for all files)
  total_lines=0

  for chr_file in ${chr_files}; do
    chr=$(basename "${chr_file}" .rsb.txt.gz | sed "s/${comp}\.//")
    echo "    Processing ${chr}..." | tee -a "${log_file}"

    # Extract data lines (skip header)
    line_count=$(zcat "${chr_file}" | tail -n +2 | tee -a "${temp_file}" | wc -l)
    total_lines=$((total_lines + line_count))

    echo "      Added ${line_count} lines" | tee -a "${log_file}"
  done

  # Compress merged file
  gzip -c "${temp_file}" > "${merged_file}"
  rm -f "${temp_file}"

  # Summary for this comparison
  echo "" | tee -a "${log_file}"
  echo "  Created: ${merged_file}" | tee -a "${log_file}"
  echo "  Total data lines: ${total_lines}" | tee -a "${log_file}"
  echo "  File size: $(du -h ${merged_file} | cut -f1)" | tee -a "${log_file}"

  # Show first few lines
  echo "" | tee -a "${log_file}"
  echo "  First 5 data lines:" | tee -a "${log_file}"
  zcat "${merged_file}" | head -n 6 | tail -n 5 | tee -a "${log_file}"

  echo "" | tee -a "${log_file}"
done

# ========================================
# Overall Summary
# ========================================

echo "========================================" | tee -a "${log_file}"
echo "SUMMARY" | tee -a "${log_file}"
echo "========================================" | tee -a "${log_file}"

ines_count=$(ls -1 ${outdir}/*.genome_wide_ines.txt.gz 2>/dev/null | wc -l)
rsb_count=$(ls -1 ${outdir}/*.genome_wide_rsb.txt.gz 2>/dev/null | wc -l)

echo "inES/iHS files created: ${ines_count} / ${#POPS[@]}" | tee -a "${log_file}"
echo "Rsb files created: ${rsb_count} / ${#COMPARISONS[@]}" | tee -a "${log_file}"
echo "" | tee -a "${log_file}"

total_expected=$((${#POPS[@]} + ${#COMPARISONS[@]}))
total_created=$((${ines_count} + ${rsb_count}))

if [[ ${total_created} -eq ${total_expected} ]]; then
  echo "SUCCESS: All genome-wide files created!" | tee -a "${log_file}"
  echo "" | tee -a "${log_file}"
  echo "Merged files:" | tee -a "${log_file}"
  ls -lh ${outdir}/*.genome_wide_*.txt.gz | tee -a "${log_file}"
  echo "" | tee -a "${log_file}"
  echo "Next step: Rscript code/12_Sweep_Detection/ehh/04_extract_peak_ehh.R" | tee -a "${log_file}"
else
  echo "WARNING: Not all merged files were created (${total_created} / ${total_expected})" | tee -a "${log_file}"
  echo "Check that all scan jobs completed successfully" | tee -a "${log_file}"
fi

echo "" | tee -a "${log_file}"
echo "Log file: ${log_file}" | tee -a "${log_file}"
