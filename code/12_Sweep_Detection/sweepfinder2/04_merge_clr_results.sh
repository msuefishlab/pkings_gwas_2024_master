#!/bin/bash
set -euo pipefail

## 04_merge_clr_results.sh
## Merge per-chromosome CLR results into genome-wide files
##
## This script combines SweepFinder2 output from all chromosomes
## for each population into single genome-wide CLR files.
##
## Usage:
##   bash code/12_Sweep_Detection/04_merge_clr_results.sh
##
## Prerequisites:
##   - 03_submit_sweepfinder2.sh jobs completed
##
## Author: Jason Gallant Lab
## Date: 2024

root="$(git rev-parse --show-toplevel)"
source "${root}/pkings_gwas.env"

# Input and output directories
indir="${root}/output_data/12_Sweep_Detection/sweepfinder2/clr_results"
outdir="${root}/output_data/12_Sweep_Detection/sweepfinder2/merged"
mkdir -p "${outdir}"

# Check if input directory exists
if [[ ! -d "${indir}" ]]; then
  echo "ERROR: CLR results directory not found: ${indir}"
  echo "Please run 03_submit_sweepfinder2.sh and wait for jobs to complete"
  exit 1
fi

# Population sets
POPS=(BP1 TP1 BP2 TP2 BP3)

# Log file
log_file="${outdir}/merge.log"
echo "Starting CLR results merge: $(date)" > "${log_file}"

echo ""
echo "=========================================="
echo "Merging CLR results"
echo "=========================================="
echo ""

# Process each population
for pop in "${POPS[@]}"; do
  echo "========================================" | tee -a "${log_file}"
  echo "Merging CLR results for ${pop}" | tee -a "${log_file}"
  echo "========================================" | tee -a "${log_file}"

  # Find all chromosome CLR files for this population
  clr_files=(${indir}/${pop}.chr*.clr.txt)

  if [[ ! -e "${clr_files[0]}" ]]; then
    echo "WARNING: No CLR files found for ${pop}" | tee -a "${log_file}"
    echo "  Pattern: ${indir}/${pop}.chr*.clr.txt" | tee -a "${log_file}"
    continue
  fi

  echo "Found ${#clr_files[@]} chromosome CLR files" | tee -a "${log_file}"

  # Output merged file
  merged_file="${outdir}/${pop}.genome_wide_clr.txt"

  # SweepFinder2 output format:
  # location    LR      alpha
  # 5000        2.34    0.002
  #
  # We want to add chromosome column:
  # chromosome  position  CLR  alpha

  echo "Creating merged file: ${merged_file}" | tee -a "${log_file}"

  # Write header
  echo -e "chromosome\tposition\tCLR\talpha" > "${merged_file}"

  # Merge all chromosomes
  for chr_file in "${clr_files[@]}"; do
    # Extract chromosome from filename
    # Example: BP1.chr6.clr.txt -> chr6
    chr=$(basename "${chr_file}" .clr.txt | sed "s/${pop}\.//")

    echo "  Merging ${chr}" | tee -a "${log_file}"

    # Append chromosome data (skip header, add chromosome column)
    # SweepFinder2 format: location LR alpha
    # Output format: chromosome location LR alpha
    tail -n +2 "${chr_file}" | awk -v chr="${chr}" 'BEGIN{OFS="\t"}{print chr, $1, $2, $3}' \
      >> "${merged_file}"
  done

  # Summary statistics
  if [[ -f "${merged_file}" ]]; then
    total_lines=$(tail -n +2 "${merged_file}" | wc -l)
    file_size=$(du -h "${merged_file}" | cut -f1)
    echo ""  | tee -a "${log_file}"
    echo "Merged file created: ${merged_file}" | tee -a "${log_file}"
    echo "  Total positions: ${total_lines}" | tee -a "${log_file}"
    echo "  File size: ${file_size}" | tee -a "${log_file}"
    echo ""  | tee -a "${log_file}"

    # Calculate CLR statistics
    echo "CLR statistics:" | tee -a "${log_file}"
    awk 'NR>1{sum+=$3; sumsq+=$3*$3; if($3>max)max=$3}
         END{
           mean=sum/NR;
           sd=sqrt(sumsq/NR - mean*mean);
           printf "  Mean CLR: %.2f\n  SD: %.2f\n  Max CLR: %.2f\n", mean, sd, max
         }' "${merged_file}" | tee -a "${log_file}"
  else
    echo "ERROR: Merged file not created!" | tee -a "${log_file}"
  fi

  echo "" | tee -a "${log_file}"
done

echo ""
echo "========================================" | tee -a "${log_file}"
echo "Merge complete: $(date)" | tee -a "${log_file}"
echo "========================================" | tee -a "${log_file}"
echo "" | tee -a "${log_file}"
echo "Output directory: ${outdir}/" | tee -a "${log_file}"
echo "Log file: ${log_file}" | tee -a "${log_file}"
echo "" | tee -a "${log_file}"

# Summary
echo "Merged files created:" | tee -a "${log_file}"
for pop in "${POPS[@]}"; do
  merged_file="${outdir}/${pop}.genome_wide_clr.txt"
  if [[ -f "${merged_file}" ]]; then
    lines=$(tail -n +2 "${merged_file}" | wc -l)
    size=$(du -h "${merged_file}" | cut -f1)
    echo "  ${pop}: ${lines} positions, ${size}" | tee -a "${log_file}"
  else
    echo "  ${pop}: NOT CREATED" | tee -a "${log_file}"
  fi
done

echo "" | tee -a "${log_file}"
echo "Next step: Run 05_extract_peak_clr.R to extract GWAS peak CLR values" | tee -a "${log_file}"
