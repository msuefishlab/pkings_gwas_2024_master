#!/bin/bash
set -euo pipefail

## 01_submit_ines_scans.sh
## Submit SLURM jobs for genome-wide inES/iHS scans
## One job per population × chromosome combination
##
## Usage:
##   bash code/12_Sweep_Detection/ehh/01_submit_ines_scans.sh [--all-chromosomes]
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
ines_dir="${outdir}/ines_scans"
log_dir="${root}/output_data/slurm_logs/12_Sweep_Detection/ehh"
mkdir -p "${ines_dir}"
mkdir -p "${log_dir}"

# Population sets to process
POPS=(BP1 BP2 BP3 TP1 TP2)

# Chromosomes to process
if [[ "${ALL_CHROM}" == "true" ]]; then
  echo "Processing all 25 chromosomes"
  CHROMOSOMES=(chr{1..25})
else
  echo "Processing GWAS peak chromosomes only (chr6, chr8, chr13, chr16, chr17, chr24)"
  CHROMOSOMES=(chr6 chr8 chr13 chr16 chr17 chr24)
fi

# Calculate total jobs
total_jobs=$((${#POPS[@]} * ${#CHROMOSOMES[@]}))

echo "========================================="
echo "inES/iHS Genome-Wide Scan Job Submission"
echo "========================================="
echo "Populations: ${#POPS[@]} (${POPS[@]})"
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
for pop in "${POPS[@]}"; do
  echo ""
  echo "Population: ${pop}"
  echo "-----------------------------------"

  for chr in "${CHROMOSOMES[@]}"; do
    # Check if polarized VCF exists
    vcf_full="${root}/output_data/12_Sweep_Detection/${pop}.polarized.vcf.gz"

    if [[ ! -f "${vcf_full}" ]]; then
      echo "  [SKIP] ${chr}: No polarized VCF found (${vcf_full})"
      skipped=$((skipped + 1))
      continue
    fi

    # Job name
    job_name="ines_${pop}_${chr}"

    # Log file
    log_file="${log_dir}/${job_name}.%j.log"

    # Submit SLURM job
    sbatch --job-name="${job_name}" \
           --output="${log_file}" \
           --export=ALL,root=${root},POP=${pop},CHR=${chr},VCF=${vcf_full} \
           "${root}/code/12_Sweep_Detection/ehh/01_run_ines_scan.sb"

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
echo "Cancel all jobs with: scancel -u \$USER -n ines"
echo ""
echo "Logs will be written to: ${log_dir}/ines_*.log"
echo ""
echo "After jobs complete, run: bash code/12_Sweep_Detection/ehh/02_submit_rsb_scans.sh"
echo "========================================="
