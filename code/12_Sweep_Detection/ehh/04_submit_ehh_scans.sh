#!/bin/bash
set -euo pipefail

## 04_submit_ehh_scans.sh
## Submit SLURM jobs for genome-wide integrated EHH scans (paired populations)
## One job per pair × chromosome combination
##
## This script replaces:
##   - DEPRECATED_submit_ines_scans.sh (deprecated)
##   - DEPRECATED_submit_rsb_scans.sh (deprecated)
##
## Usage:
##   bash code/12_Sweep_Detection/ehh/04_submit_ehh_scans.sh [--all-chromosomes]
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
outdir="${root}/output_data/12_Sweep_Detection/ehh"
scan_dir="${outdir}/scans"
log_dir="${root}/output_data/slurm_logs/12_Sweep_Detection/ehh"
mkdir -p "${scan_dir}"
mkdir -p "${log_dir}"

# ========================================
# Define paired comparisons
# ========================================

# Each pair has: pair name, BP list file, TP list file
declare -A PAIR_BP_LIST
declare -A PAIR_TP_LIST

PAIRS=(BP1_TP1 BP2_TP2 BP3_TP2)

PAIR_BP_LIST[BP1_TP1]="${root}/input_data/12_Sweep_Detection/split_up/BP1.txt"
PAIR_TP_LIST[BP1_TP1]="${root}/input_data/12_Sweep_Detection/split_up/TP1.txt"

PAIR_BP_LIST[BP2_TP2]="${root}/input_data/12_Sweep_Detection/split_up/BP2.txt"
PAIR_TP_LIST[BP2_TP2]="${root}/input_data/12_Sweep_Detection/split_up/TP2.txt"

PAIR_BP_LIST[BP3_TP2]="${root}/input_data/12_Sweep_Detection/split_up/BP3.txt"
PAIR_TP_LIST[BP3_TP2]="${root}/input_data/12_Sweep_Detection/split_up/TP2.txt"

# ========================================
# Chromosomes to process
# ========================================

if [[ "${ALL_CHROM}" == "true" ]]; then
  echo "Processing all 25 chromosomes"
  CHROMOSOMES=(chr{1..25})
else
  echo "Processing GWAS peak chromosomes only (chr6, chr8, chr13, chr16, chr17, chr24)"
  CHROMOSOMES=(chr6 chr8 chr13 chr16 chr17 chr24)
fi

# ========================================
# Calculate total jobs
# ========================================

total_jobs=$((${#PAIRS[@]} * ${#CHROMOSOMES[@]}))

echo "========================================="
echo "Integrated EHH Scan Job Submission"
echo "========================================="
echo "Pairs:       ${#PAIRS[@]} (${PAIRS[@]})"
echo "Chromosomes: ${#CHROMOSOMES[@]} (${CHROMOSOMES[@]})"
echo "Total jobs:  ${total_jobs}"
echo "========================================="
echo ""

# ========================================
# Validate prerequisites
# ========================================

echo "Checking prerequisites..."
echo ""

vcf_dir="${root}/output_data/12_Sweep_Detection/by_chrom"
all_vcfs_exist=true

for pair in "${PAIRS[@]}"; do
  # Check one chromosome as a test
  test_vcf="${vcf_dir}/${pair}.polarized.chr6.vcf.gz"
  if [[ ! -f "${test_vcf}" ]]; then
    echo "  [WARNING] Missing VCF for ${pair} (tested ${test_vcf})"
    echo "            Run 03_split_by_chromosome.sh first!"
    all_vcfs_exist=false
  else
    echo "  [OK] ${pair} VCFs exist"
  fi

  # Check population lists
  bp_list="${PAIR_BP_LIST[${pair}]}"
  tp_list="${PAIR_TP_LIST[${pair}]}"

  if [[ ! -f "${bp_list}" ]]; then
    echo "  [ERROR] BP list not found: ${bp_list}"
    all_vcfs_exist=false
  fi

  if [[ ! -f "${tp_list}" ]]; then
    echo "  [ERROR] TP list not found: ${tp_list}"
    all_vcfs_exist=false
  fi
done

echo ""

if [[ "${all_vcfs_exist}" == "false" ]]; then
  echo "ERROR: Missing required files. Cannot proceed."
  echo ""
  echo "To split VCFs by chromosome, run:"
  echo "  bash code/12_Sweep_Detection/ehh/03_split_by_chromosome.sh"
  exit 1
fi

# ========================================
# Confirm submission
# ========================================

read -p "Submit ${total_jobs} jobs? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
  echo "Cancelled."
  exit 0
fi

# ========================================
# Submit jobs
# ========================================

submitted=0
skipped=0

for pair in "${PAIRS[@]}"; do
  echo ""
  echo "Pair: ${pair}"
  echo "-----------------------------------"

  bp_list="${PAIR_BP_LIST[${pair}]}"
  tp_list="${PAIR_TP_LIST[${pair}]}"

  for chr in "${CHROMOSOMES[@]}"; do
    # VCF file for this pair and chromosome
    vcf_file="${vcf_dir}/${pair}.polarized.${chr}.vcf.gz"

    if [[ ! -f "${vcf_file}" ]]; then
      echo "  [SKIP] ${chr}: VCF not found (${vcf_file})"
      skipped=$((skipped + 1))
      continue
    fi

    # Job name
    job_name="ehh_${pair}_${chr}"

    # Log file
    log_file="${log_dir}/${job_name}.%j.log"

    # Submit SLURM job
    sbatch --job-name="${job_name}" \
           --output="${log_file}" \
           --export=ALL,root=${root},PAIR=${pair},CHR=${chr},VCF=${vcf_file},BP_LIST=${bp_list},TP_LIST=${tp_list},OUTPUT_DIR=${scan_dir} \
           "${root}/code/12_Sweep_Detection/ehh/04_run_ehh_scan.sb"

    echo "  [SUBMITTED] ${chr}: ${job_name}"
    submitted=$((submitted + 1))
  done
done

echo ""
echo "========================================="
echo "SUBMISSION SUMMARY"
echo "========================================="
echo "Jobs submitted: ${submitted}"
echo "Jobs skipped:   ${skipped}"
echo "Total:          ${total_jobs}"
echo ""
echo "Monitor jobs with: squeue -u \$USER"
echo "Cancel all jobs with: scancel -u \$USER -n ehh"
echo ""
echo "Logs will be written to: ${log_dir}/ehh_*.log"
echo ""
echo "After jobs complete, run:"
echo "  bash code/12_Sweep_Detection/ehh/05_merge_ehh_results.sh"
echo "========================================="
