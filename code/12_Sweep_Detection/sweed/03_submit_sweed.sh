#!/bin/bash
set -euo pipefail

## 03_submit_sweed.sh
## Submit SLURM jobs for SweeD CLR scan
## One job per population × chromosome combination
##
## Usage:
##   bash code/12_Sweep_Detection/sweed/03_submit_sweed.sh [--all-chromosomes]
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
outdir="${root}/output_data/12_Sweep_Detection/sweed"
clr_dir="${outdir}/clr_results"
log_dir="${root}/output_data/slurm_logs/12_Sweep_Detection"
mkdir -p "${clr_dir}"
mkdir -p "${log_dir}"

# Check if SweeD container exists
if [[ ! -f "${sweed_image}" ]]; then
  echo "ERROR: SweeD container not found: ${sweed_image}"
  echo "Please build the container first."
  exit 1
fi

# Population sets to process
POPS=(BP1 TP1 BP2 TP2 BP3)

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
echo "SweeD CLR Scan Job Submission"
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
    # Check if VCF exists (per-chromosome or full)
    vcf_by_chr="${root}/output_data/12_Sweep_Detection/vcfs/by_chrom/${pop}.polarized.${chr}.vcf.gz"
    vcf_full="${root}/output_data/12_Sweep_Detection/vcfs/${pop}.polarized.vcf.gz"

    if [[ ! -f "${vcf_by_chr}" ]] && [[ ! -f "${vcf_full}" ]]; then
      echo "  [SKIP] ${chr}: No VCF found"
      skipped=$((skipped + 1))
      continue
    fi

    # Set VCF path for job
    if [[ -f "${vcf_by_chr}" ]]; then
      vcf="${vcf_by_chr}"
    else
      vcf="${vcf_full}"
    fi

    # Job name
    job_name="sweed_${pop}_${chr}"

    # Log file
    log_file="${log_dir}/${job_name}.%j.log"

    # Submit SLURM job
    sbatch --job-name="${job_name}" \
           --output="${log_file}" \
           --export=ALL,root=${root},POP=${pop},CHR=${chr},VCF=${vcf} \
           "${root}/code/12_Sweep_Detection/sweed/03_run_sweed.sb"

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
echo "Cancel all jobs with: scancel -u \$USER -n sweed"
echo ""
echo "Logs will be written to: ${log_dir}/sweed_*.log"
echo ""
echo "After jobs complete, run: code/12_Sweep_Detection/sweed/04_merge_clr_results.sh"
echo "========================================="
