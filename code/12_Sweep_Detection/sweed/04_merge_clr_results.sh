#!/bin/bash
set -euo pipefail

## 04_merge_clr_results.sh
## Merge per-chromosome SweeD CLR results into genome-wide files
##
## SweeD output format (with // comment lines):
##   // Position    Likelihood  Alpha
##   1000           2.456       0.0015
##   2000           1.234       0.0020
##
## Target merged format (for R script compatibility):
##   chromosome  position  CLR  alpha
##   chr6        1000      2.456  0.0015
##   chr6        2000      1.234  0.0020
##
## Usage:
##   bash code/12_Sweep_Detection/sweed/04_merge_clr_results.sh
##
## Author: Jason Gallant Lab
## Date: 2024

root="$(git rev-parse --show-toplevel)"
source "${root}/pkings_gwas.env"

# Define directories
indir="${root}/output_data/12_Sweep_Detection/sweed/clr_results"
outdir="${root}/output_data/12_Sweep_Detection/sweed/merged"
mkdir -p "${outdir}"

# Population sets
POPS=(BP1 TP1 BP2 TP2 BP3)

# Log file
log_file="${outdir}/merge.log"
echo "SweeD CLR Results Merging" > "${log_file}"
echo "Generated: $(date)" >> "${log_file}"
echo "========================================" >> "${log_file}"
echo "" >> "${log_file}"

echo "Merging per-chromosome CLR results into genome-wide files..."
echo ""

# Process each population
for pop in "${POPS[@]}"; do
  echo "========================================" | tee -a "${log_file}"
  echo "Population: ${pop}" | tee -a "${log_file}"
  echo "========================================" | tee -a "${log_file}"

  # Output merged file
  merged_file="${outdir}/${pop}.genome_wide_clr.txt"

  # Write header
  echo -e "chromosome\tposition\tCLR\talpha" > "${merged_file}"

  # Find all chromosome files for this population
  chr_files=$(ls ${indir}/${pop}.chr*.sweed.txt 2>/dev/null || true)

  if [[ -z "${chr_files}" ]]; then
    echo "  WARNING: No CLR files found for ${pop}" | tee -a "${log_file}"
    continue
  fi

  # Count files
  file_count=$(echo "${chr_files}" | wc -w)
  echo "  Found ${file_count} chromosome files" | tee -a "${log_file}"

  # Merge all chromosomes
  total_lines=0

  for chr_file in ${chr_files}; do
    # Extract chromosome from filename
    # e.g., BP1.chr6.sweed.txt -> chr6
    chr=$(basename "${chr_file}" .sweed.txt | sed "s/${pop}\.//")

    echo "    Processing ${chr}..." | tee -a "${log_file}"

    # Parse SweeD output format
    # SweeD output has:
    #   - Comment lines starting with //
    #   - Header line (after comments): "Position    Likelihood  Alpha"
    #   - Data lines: position  likelihood  alpha
    #
    # We need to:
    #   1. Skip comment lines (// ...)
    #   2. Skip header line
    #   3. Add chromosome column
    #   4. Rename Likelihood -> CLR

    line_count=$(grep -v '^//' "${chr_file}" | tail -n +2 | \
      awk -v chr="${chr}" 'BEGIN{OFS="\t"}{print chr, $1, $2, $3}' | \
      tee -a "${merged_file}" | wc -l)

    total_lines=$((total_lines + line_count))

    echo "      Added ${line_count} lines" | tee -a "${log_file}"
  done

  # Summary for this population
  echo "" | tee -a "${log_file}"
  echo "  Created: ${merged_file}" | tee -a "${log_file}"
  echo "  Total lines: ${total_lines}" | tee -a "${log_file}"
  echo "  File size: $(du -h ${merged_file} | cut -f1)" | tee -a "${log_file}"

  # Show first few lines
  echo "" | tee -a "${log_file}"
  echo "  First 5 data lines:" | tee -a "${log_file}"
  head -n 6 "${merged_file}" | tail -n 5 | tee -a "${log_file}"

  # Show last few lines
  echo "" | tee -a "${log_file}"
  echo "  Last 5 data lines:" | tee -a "${log_file}"
  tail -n 5 "${merged_file}" | tee -a "${log_file}"

  echo "" | tee -a "${log_file}"
done

# Overall summary
echo "========================================" | tee -a "${log_file}"
echo "SUMMARY" | tee -a "${log_file}"
echo "========================================" | tee -a "${log_file}"

merged_count=$(ls -1 ${outdir}/*.genome_wide_clr.txt 2>/dev/null | wc -l)
echo "Merged files created: ${merged_count} / ${#POPS[@]}" | tee -a "${log_file}"
echo "" | tee -a "${log_file}"

if [[ ${merged_count} -eq ${#POPS[@]} ]]; then
  echo "SUCCESS: All genome-wide CLR files created!" | tee -a "${log_file}"
  echo "" | tee -a "${log_file}"
  echo "Merged files:" | tee -a "${log_file}"
  ls -lh ${outdir}/*.genome_wide_clr.txt | tee -a "${log_file}"
  echo "" | tee -a "${log_file}"
  echo "Next step: Rscript code/12_Sweep_Detection/sweed/05_extract_peak_clr.R" | tee -a "${log_file}"
else
  echo "WARNING: Not all merged files were created" | tee -a "${log_file}"
  echo "Check that SweeD jobs completed successfully" | tee -a "${log_file}"
fi

echo "" | tee -a "${log_file}"
echo "Log file: ${log_file}" | tee -a "${log_file}"
