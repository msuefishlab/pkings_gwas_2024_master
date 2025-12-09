#!/bin/bash
set -euo pipefail

## 02_submit_rsb_scans.sh
## Submit SLURM jobs for genome-wide Rsb scans
## One job per comparison × chromosome combination
##
## Usage:
##   bash code/12_Sweep_Detection/ehh/02_submit_rsb_scans.sh [--all-chromosomes]
##
## Prerequisites:
##   - 01_submit_ines_scans.sh must be completed (inES files must exist)
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
rsb_dir="${outdir}/rsb_scans"
log_dir="${root}/output_data/slurm_logs/12_Sweep_Detection/ehh"
mkdir -p "${rsb_dir}"
mkdir -p "${log_dir}"

# Population comparisons (BP vs TP)
# Format: "POP1 POP2"
COMPARISONS=(
  "BP1 TP1"
  "BP2 TP2"
  "BP3 TP2"
)

# Chromosomes to process
if [[ "${ALL_CHROM}" == "true" ]]; then
  echo "Processing all 25 chromosomes"
  CHROMOSOMES=(chr{1..25})
else
  echo "Processing GWAS peak chromosomes only (chr6, chr8, chr13, chr16, chr17, chr24)"
  CHROMOSOMES=(chr6 chr8 chr13 chr16 chr17 chr24)
fi

# Calculate total jobs
total_jobs=$((${#COMPARISONS[@]} * ${#CHROMOSOMES[@]}))

echo "========================================="
echo "Rsb Genome-Wide Scan Job Submission"
echo "========================================="
echo "Comparisons: ${#COMPARISONS[@]} (BP1-TP1, BP2-TP2, BP3-TP2)"
echo "Chromosomes: ${#CHROMOSOMES[@]} (${CHROMOSOMES[@]})"
echo "Total jobs: ${total_jobs}"
echo "========================================="
echo ""

# Confirm submission
read -p "Submit ${total_jobs} jobs? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
  echo "Cancelled."
  exit 0
fi

# Counter for submitted jobs
submitted=0
skipped=0

# Submit jobs
for comp in "${COMPARISONS[@]}"; do
  # Split comparison into pop1 and pop2
  read -r pop1 pop2 <<< "$comp"

  echo ""
  echo "Comparison: ${pop1} vs ${pop2}"
  echo "-----------------------------------"

  for chr in "${CHROMOSOMES[@]}"; do
    # Check if required inES files exist
    ines1="${outdir}/ines_scans/${pop1}.${chr}.ines.txt.gz"
    ines2="${outdir}/ines_scans/${pop2}.${chr}.ines.txt.gz"

    if [[ ! -f "${ines1}" ]]; then
      echo "  [SKIP] ${chr}: inES file not found for ${pop1} (${ines1})"
      skipped=$((skipped + 1))
      continue
    fi

    if [[ ! -f "${ines2}" ]]; then
      echo "  [SKIP] ${chr}: inES file not found for ${pop2} (${ines2})"
      skipped=$((skipped + 1))
      continue
    fi

    # Job name
    job_name="rsb_${pop1}_${pop2}_${chr}"

    # Log file
    log_file="${log_dir}/${job_name}.%j.log"

    # Submit SLURM job
    sbatch --job-name="${job_name}" \
           --output="${log_file}" \
           --export=ALL,root=${root},POP1=${pop1},POP2=${pop2},CHR=${chr} \
           "${root}/code/12_Sweep_Detection/ehh/02_run_rsb_scan.sb"

    echo "  [SUBMITTED] ${chr}: ${job_name}"
    submitted=$((submitted + 1))
  done
done

echo ""
echo "========================================="
echo "SUBMISSION SUMMARY"
echo "========================================="
echo "Jobs submitted: ${submitted}"
echo "Jobs skipped: ${skipped}"
echo "Total: ${total_jobs}"
echo ""
echo "Monitor jobs with: squeue -u \$USER"
echo "Cancel all jobs with: scancel -u \$USER -n rsb"
echo ""
echo "Logs will be written to: ${log_dir}/rsb_*.log"
echo ""
echo "After jobs complete, run: bash code/12_Sweep_Detection/ehh/03_merge_ehh_results.sh"
echo "========================================="
