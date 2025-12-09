#!/bin/bash
set -euo pipefail

## 03_submit_sweepfinder2.sh
## Submit SweepFinder2 jobs for all populations and chromosomes
##
## This script submits SLURM jobs to run genome-wide CLR scans.
## Each job runs SweepFinder2 for one population, one chromosome.
##
## Usage:
##   bash code/12_Sweep_Detection/sweepfinder2/03_submit_sweepfinder2.sh [--all-chromosomes]
##
## Options:
##   --all-chromosomes    Submit jobs for all 25 chromosomes instead of just GWAS peaks
##
## Prerequisites:
##   - 01_convert_vcf_to_sf2.sh completed
##   - 02_calculate_background_sfs.sh completed
##   - SweepFinder2 container built
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

# Create log directory
mkdir -p "${root}/output_data/slurm_logs/12_Sweep_Detection"

# Population sets
POPS=(BP1 TP1 BP2 TP2 BP3)

# Chromosomes to process
if [[ "${ALL_CHROM}" == "true" ]]; then
  echo "Submitting jobs for all 25 chromosomes"
  CHROMOSOMES=(chr{1..25})
else
  echo "Submitting jobs for GWAS peak chromosomes only"
  CHROMOSOMES=(chr6 chr8 chr13 chr16 chr17 chr24)
fi

# Check prerequisites
sfs_dir="${root}/output_data/12_Sweep_Detection/sweepfinder2/sfs_input"
spec_dir="${root}/output_data/12_Sweep_Detection/sweepfinder2/background"

if [[ ! -d "${sfs_dir}" ]] || [[ $(ls ${sfs_dir}/*.sf2.txt 2>/dev/null | wc -l) -eq 0 ]]; then
  echo "ERROR: SFS input files not found in ${sfs_dir}"
  echo "Please run 01_convert_vcf_to_sf2.sh first"
  exit 1
fi

if [[ ! -d "${spec_dir}" ]] || [[ $(ls ${spec_dir}/*.spectrum 2>/dev/null | wc -l) -eq 0 ]]; then
  echo "ERROR: Background spectrum files not found in ${spec_dir}"
  echo "Please run 02_calculate_background_sfs.sh first"
  exit 1
fi

# Calculate total jobs
total_jobs=$((${#POPS[@]} * ${#CHROMOSOMES[@]}))
echo ""
echo "=========================================="
echo "Submitting SweepFinder2 CLR scan jobs"
echo "=========================================="
echo "Populations: ${#POPS[@]} (${POPS[@]})"
echo "Chromosomes: ${#CHROMOSOMES[@]}"
echo "Total jobs: ${total_jobs}"
echo "=========================================="
echo ""

# Track job submissions
submitted=0
skipped=0

# Job array for tracking
declare -a job_ids

# Submit jobs
for pop in "${POPS[@]}"; do
  for chr in "${CHROMOSOMES[@]}"; do

    # Check if input files exist
    sfs_file="${sfs_dir}/${pop}.${chr}.sf2.txt"
    spec_file="${spec_dir}/${pop}.spectrum"

    if [[ ! -f "${sfs_file}" ]]; then
      echo "WARNING: Skipping ${pop} ${chr} - SFS file not found: ${sfs_file}"
      skipped=$((skipped + 1))
      continue
    fi

    if [[ ! -f "${spec_file}" ]]; then
      echo "WARNING: Skipping ${pop} ${chr} - spectrum file not found: ${spec_file}"
      skipped=$((skipped + 1))
      continue
    fi

    # Job name
    job_name="sf2_${pop}_${chr}"

    # Log file
    log_file="${root}/output_data/slurm_logs/12_Sweep_Detection/${job_name}.log"

    # Submit job
    echo "Submitting: ${pop} ${chr}"

    job_id=$(sbatch --job-name="${job_name}" \
           --output="${log_file}" \
           --export=ALL,root=${root},POP=${pop},CHR=${chr} \
           "${root}/code/12_Sweep_Detection/sweepfinder2/03_run_sweepfinder2.sb" \
           | awk '{print $NF}')

    echo "  Job ID: ${job_id}"
    echo "  Log: ${log_file}"

    job_ids+=("${job_id}")
    submitted=$((submitted + 1))
  done
done

echo ""
echo "=========================================="
echo "Job submission complete"
echo "=========================================="
echo "Jobs submitted: ${submitted}"
echo "Jobs skipped: ${skipped}"
echo "Total: $((submitted + skipped)) / ${total_jobs}"
echo ""

if [[ ${submitted} -gt 0 ]]; then
  echo "Monitor jobs with:"
  echo "  squeue -u \$USER"
  echo "  squeue -j $(IFS=,; echo "${job_ids[*]}")"
  echo ""
  echo "View log files in:"
  echo "  ${root}/output_data/slurm_logs/12_Sweep_Detection/"
  echo ""
  echo "After jobs complete, run:"
  echo "  bash code/12_Sweep_Detection/sweepfinder2/04_merge_clr_results.sh"
fi

# Save job IDs to file for later reference
job_ids_file="${root}/output_data/12_Sweep_Detection/sweepfinder2/submitted_job_ids.txt"
printf "%s\n" "${job_ids[@]}" > "${job_ids_file}"
echo "Job IDs saved to: ${job_ids_file}"
