#!/bin/bash
set -euo pipefail

## 05_merge_ehh_results.sh
## Merge per-chromosome inES/iHS and Rsb results into genome-wide files
##
## Usage:
##   bash code/12_Sweep_Detection/ehh/05_merge_ehh_results.sh
##
## Prerequisites:
##   - 04_submit_ehh_scans.sh must be completed (integrated scan files must exist)
##
## Author: Jason Gallant Lab
## Date: 2024

root="$(git rev-parse --show-toplevel)"
source "${root}/pkings_gwas.env"

# Define directories
scan_dir="${root}/output_data/12_Sweep_Detection/ehh/scans"
outdir="${root}/output_data/12_Sweep_Detection/ehh/merged"
mkdir -p "${outdir}"

# Paired populations
PAIRS=(BP1_TP1 BP2_TP2 BP3_TP2)

# Log file
log_file="${outdir}/merge.log"
echo "EHH Results Merging (Paired Population Pipeline)" > "${log_file}"
echo "Generated: $(date)" >> "${log_file}"
echo "========================================" >> "${log_file}"
echo "" >> "${log_file}"

echo "========================================"
echo "Merging per-chromosome EHH results"
echo "========================================"
echo ""

# ========================================
# Part 1: Merge inES/iHS files (BP populations)
# ========================================

echo "Part 1a: Merging BP population inES/iHS results..." | tee -a "${log_file}"
echo "" | tee -a "${log_file}"

for pair in "${PAIRS[@]}"; do
  # Extract BP label from pair (e.g., BP1_TP1 -> BP1)
  bp_label=$(echo "${pair}" | cut -d'_' -f1)

  echo "========================================" | tee -a "${log_file}"
  echo "Pair: ${pair} (BP population: ${bp_label})" | tee -a "${log_file}"
  echo "========================================" | tee -a "${log_file}"

  # Output merged file
  merged_file="${outdir}/${bp_label}.genome_wide_ines.txt.gz"

  # Find all chromosome files for this pair's BP population
  chr_files=$(ls ${scan_dir}/${pair}.BP.chr*.ines.txt.gz 2>/dev/null || true)

  if [[ -z "${chr_files}" ]]; then
    echo "  WARNING: No BP inES files found for ${pair}" | tee -a "${log_file}"
    echo "" | tee -a "${log_file}"
    continue
  fi

  # Count files
  file_count=$(echo "${chr_files}" | wc -w)
  echo "  Found ${file_count} chromosome files" | tee -a "${log_file}"

  # Create temporary file for merging
  temp_file="${outdir}/.temp_${bp_label}_ines.txt"

  # Process first file to get header
  first_file=$(echo "${chr_files}" | awk '{print $1}')
  zcat "${first_file}" | head -n 1 > "${temp_file}"

  # Merge all chromosomes (skip headers for all files)
  total_lines=0

  for chr_file in ${chr_files}; do
    chr=$(basename "${chr_file}" .ines.txt.gz | sed "s/${pair}\.BP\.//")
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
# Part 1b: Merge inES/iHS files (TP populations)
# ========================================

echo "" | tee -a "${log_file}"
echo "Part 1b: Merging TP population inES/iHS results..." | tee -a "${log_file}"
echo "" | tee -a "${log_file}"

for pair in "${PAIRS[@]}"; do
  # Extract TP label from pair (e.g., BP1_TP1 -> TP1)
  tp_label=$(echo "${pair}" | cut -d'_' -f2)

  echo "========================================" | tee -a "${log_file}"
  echo "Pair: ${pair} (TP population: ${tp_label})" | tee -a "${log_file}"
  echo "========================================" | tee -a "${log_file}"

  # Output merged file
  merged_file="${outdir}/${tp_label}.genome_wide_ines.txt.gz"

  # Find all chromosome files for this pair's TP population
  chr_files=$(ls ${scan_dir}/${pair}.TP.chr*.ines.txt.gz 2>/dev/null || true)

  if [[ -z "${chr_files}" ]]; then
    echo "  WARNING: No TP inES files found for ${pair}" | tee -a "${log_file}"
    echo "" | tee -a "${log_file}"
    continue
  fi

  # Count files
  file_count=$(echo "${chr_files}" | wc -w)
  echo "  Found ${file_count} chromosome files" | tee -a "${log_file}"

  # Create temporary file for merging
  temp_file="${outdir}/.temp_${tp_label}_ines.txt"

  # Check if this is first time processing this TP population
  # (TP2 appears in both BP2_TP2 and BP3_TP2)
  if [[ -f "${merged_file}" ]]; then
    echo "  NOTE: ${tp_label} already processed (appears in multiple pairs)" | tee -a "${log_file}"
    echo "        Deduplicating positions..." | tee -a "${log_file}"

    # Append to existing file (will deduplicate later)
    existing_merged=true
    zcat "${merged_file}" > "${temp_file}"
  else
    existing_merged=false
    # Process first file to get header
    first_file=$(echo "${chr_files}" | awk '{print $1}')
    zcat "${first_file}" | head -n 1 > "${temp_file}"
  fi

  # Merge all chromosomes (skip headers for all files)
  total_lines=0

  for chr_file in ${chr_files}; do
    chr=$(basename "${chr_file}" .ines.txt.gz | sed "s/${pair}\.TP\.//")
    echo "    Processing ${chr}..." | tee -a "${log_file}"

    # Extract data lines (skip header)
    line_count=$(zcat "${chr_file}" | tail -n +2 | tee -a "${temp_file}" | wc -l)
    total_lines=$((total_lines + line_count))

    echo "      Added ${line_count} lines" | tee -a "${log_file}"
  done

  # If TP population appears in multiple pairs, deduplicate by CHR:POSITION
  if [[ "${tp_label}" == "TP2" ]] || [[ "${existing_merged}" == "true" ]]; then
    echo "  Deduplicating positions (keeping first occurrence)..." | tee -a "${log_file}"
    temp_dedup="${outdir}/.temp_${tp_label}_ines_dedup.txt"
    awk 'NR==1 || !seen[$1":"$2]++ {print}' "${temp_file}" > "${temp_dedup}"
    mv "${temp_dedup}" "${temp_file}"

    # Count final unique positions
    final_lines=$(($(wc -l < "${temp_file}") - 1))  # Subtract header
    echo "  After deduplication: ${final_lines} unique positions" | tee -a "${log_file}"
  fi

  # Compress merged file
  gzip -c "${temp_file}" > "${merged_file}"
  rm -f "${temp_file}"

  # Summary for this population
  echo "" | tee -a "${log_file}"
  echo "  Created: ${merged_file}" | tee -a "${log_file}"
  echo "  Total data lines: $(zcat ${merged_file} | tail -n +2 | wc -l)" | tee -a "${log_file}"
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

for pair in "${PAIRS[@]}"; do
  echo "========================================" | tee -a "${log_file}"
  echo "Pair: ${pair}" | tee -a "${log_file}"
  echo "========================================" | tee -a "${log_file}"

  # Output merged file
  merged_file="${outdir}/${pair}.genome_wide_rsb.txt.gz"

  # Find all chromosome files for this pair's Rsb
  chr_files=$(ls ${scan_dir}/${pair}.chr*.rsb.txt.gz 2>/dev/null || true)

  if [[ -z "${chr_files}" ]]; then
    echo "  WARNING: No Rsb files found for ${pair}" | tee -a "${log_file}"
    echo "" | tee -a "${log_file}"
    continue
  fi

  # Count files
  file_count=$(echo "${chr_files}" | wc -w)
  echo "  Found ${file_count} chromosome files" | tee -a "${log_file}"

  # Create temporary file for merging
  temp_file="${outdir}/.temp_${pair}_rsb.txt"

  # Process first file to get header
  first_file=$(echo "${chr_files}" | awk '{print $1}')
  zcat "${first_file}" | head -n 1 > "${temp_file}"

  # Merge all chromosomes (skip headers for all files)
  total_lines=0

  for chr_file in ${chr_files}; do
    chr=$(basename "${chr_file}" .rsb.txt.gz | sed "s/${pair}\.//")
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

# Expected: 5 inES files (BP1, BP2, BP3, TP1, TP2) + 3 Rsb files (BP1_TP1, BP2_TP2, BP3_TP2)
expected_ines=5
expected_rsb=3

echo "inES/iHS files created: ${ines_count} / ${expected_ines}" | tee -a "${log_file}"
echo "Rsb files created: ${rsb_count} / ${expected_rsb}" | tee -a "${log_file}"
echo "" | tee -a "${log_file}"

total_expected=$((${expected_ines} + ${expected_rsb}))
total_created=$((${ines_count} + ${rsb_count}))

if [[ ${total_created} -eq ${total_expected} ]]; then
  echo "SUCCESS: All genome-wide files created!" | tee -a "${log_file}"
  echo "" | tee -a "${log_file}"
  echo "Merged files:" | tee -a "${log_file}"
  ls -lh ${outdir}/*.genome_wide_*.txt.gz | tee -a "${log_file}"
  echo "" | tee -a "${log_file}"
  echo "Next step: singularity exec --bind \${root}:/project_root \${rehh_image} Rscript /project_root/code/12_Sweep_Detection/ehh/04_extract_peak_ehh.R" | tee -a "${log_file}"
else
  echo "WARNING: Not all merged files were created (${total_created} / ${total_expected})" | tee -a "${log_file}"
  echo "Check that all scan jobs completed successfully" | tee -a "${log_file}"
fi

echo "" | tee -a "${log_file}"
echo "Log file: ${log_file}" | tee -a "${log_file}"
