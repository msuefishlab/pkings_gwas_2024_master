#!/bin/bash
set -euo pipefail

## 02_calculate_background_sfs.sh
## Calculate genome-wide empirical site frequency spectrum (SFS) per population
## This creates the background null model for SweeD CLR calculations
##
## Usage:
##   bash code/12_Sweep_Detection/sweed/02_calculate_background_sfs.sh
##
## Note: This uses the full polarized VCF per population to calculate genome-wide SFS
##
## Author: Jason Gallant Lab
## Date: 2024

root="$(git rev-parse --show-toplevel)"
source "${root}/pkings_gwas.env"

# Create output directories
outdir="${root}/output_data/12_Sweep_Detection/sweed"
sfs_dir="${outdir}/background"
mkdir -p "${sfs_dir}"

# Check if SweeD container exists
if [[ ! -f "${sweed_image}" ]]; then
  echo "ERROR: SweeD container not found: ${sweed_image}"
  echo "Please build the container first:"
  echo "  cd images/"
  echo "  docker build --platform linux/amd64 -f Dockerfile.sweed -t sweed ."
  echo "  docker push <username>/sweed"
  echo "  singularity build sweed.sif docker://<username>/sweed"
  exit 1
fi

# Population sets to process
POPS=(BP1 TP1 BP2 TP2 BP3)

# Log file
log_file="${sfs_dir}/sfs_calculation.log"
echo "SweeD Background SFS Calculation" > "${log_file}"
echo "Generated: $(date)" >> "${log_file}"
echo "========================================" >> "${log_file}"
echo "" >> "${log_file}"

echo "Calculating genome-wide empirical SFS for each population..."
echo "This establishes the background null model for CLR calculations."
echo ""

# Process each population
for pop in "${POPS[@]}"; do
  echo "========================================" | tee -a "${log_file}"
  echo "Population: ${pop}" | tee -a "${log_file}"
  echo "========================================" | tee -a "${log_file}"

  # Input VCF (full genome, polarized)
  vcf="${root}/output_data/12_Sweep_Detection/vcfs/by_group/${pop}.polarized.vcf.gz"

  # Check if VCF exists
  if [[ ! -f "${vcf}" ]]; then
    echo "  ERROR: VCF not found: ${vcf}" | tee -a "${log_file}"
    echo "  Run polarize_vcf.sh and split_vcfs.sh first." | tee -a "${log_file}"
    continue
  fi

  # Output SFS file
  sfs_file="${sfs_dir}/${pop}.genome_sfs.txt"

  # Temporary output directory for SweeD run
  temp_dir="${sfs_dir}/temp_${pop}"
  mkdir -p "${temp_dir}"

  echo "  Input VCF: ${vcf}" | tee -a "${log_file}"
  echo "  Output SFS: ${sfs_file}" | tee -a "${log_file}"
  echo "  Started: $(date)" | tee -a "${log_file}"
  echo "" | tee -a "${log_file}"

  # SweeD does not support gzipped VCF files directly
  # Decompress to temporary uncompressed VCF
  temp_vcf="${temp_dir}/${pop}.polarized.vcf"

  echo "  Decompressing VCF (SweeD requires uncompressed VCF)..." | tee -a "${log_file}"
  gunzip -c "${vcf}" > "${temp_vcf}"

  if [[ ! -f "${temp_vcf}" ]]; then
    echo "  ERROR: Failed to decompress VCF" | tee -a "${log_file}"
    exit_code=1
  else
    echo "  Decompressed VCF size: $(du -h ${temp_vcf} | cut -f1)" | tee -a "${log_file}"
    echo "" | tee -a "${log_file}"

    # Run SweeD with -osfs flag to output empirical SFS
    # We run on the full genome to get comprehensive background spectrum
    # -osfs: output site frequency spectrum
    # -name: run name (for output file naming)

    cd "${temp_dir}"

    singularity exec --bind ${temp_dir}:/temp,${sfs_dir}:/output ${sweed_image} \
      SweeD -input "/temp/${pop}.polarized.vcf" \
            -name "${pop}_genome_sfs" \
            -osfs "/output/${pop}.genome_sfs.txt" \
            2>&1 | tee -a "${log_file}"

    exit_code=$?

    cd "${root}"
  fi

  if [[ ${exit_code} -eq 0 ]] && [[ -f "${sfs_file}" ]]; then
    echo "" | tee -a "${log_file}"
    echo "  SUCCESS: SFS calculated for ${pop}" | tee -a "${log_file}"
    echo "  Output file: ${sfs_file}" | tee -a "${log_file}"
    echo "  File size: $(du -h ${sfs_file} | cut -f1)" | tee -a "${log_file}"
    echo "  Completed: $(date)" | tee -a "${log_file}"

    # Clean up temp directory
    rm -rf "${temp_dir}"
  else
    echo "" | tee -a "${log_file}"
    echo "  ERROR: SFS calculation failed for ${pop}" | tee -a "${log_file}"
    echo "  Exit code: ${exit_code}" | tee -a "${log_file}"
    echo "  Check logs above for details" | tee -a "${log_file}"
  fi

  echo "" | tee -a "${log_file}"
done

# Summary
echo "========================================" | tee -a "${log_file}"
echo "SUMMARY" | tee -a "${log_file}"
echo "========================================" | tee -a "${log_file}"

sfs_count=$(ls -1 ${sfs_dir}/*.genome_sfs.txt 2>/dev/null | wc -l)
echo "SFS files created: ${sfs_count} / ${#POPS[@]}" | tee -a "${log_file}"
echo "" | tee -a "${log_file}"

if [[ ${sfs_count} -eq ${#POPS[@]} ]]; then
  echo "SUCCESS: All background SFS files calculated!" | tee -a "${log_file}"
  echo "" | tee -a "${log_file}"
  echo "SFS files:" | tee -a "${log_file}"
  ls -lh ${sfs_dir}/*.genome_sfs.txt | tee -a "${log_file}"
  echo "" | tee -a "${log_file}"
  echo "Next step: Run 03_submit_sweed.sh to calculate CLR statistics" | tee -a "${log_file}"
else
  echo "WARNING: Not all SFS files were created successfully" | tee -a "${log_file}"
  echo "Check the log above for errors" | tee -a "${log_file}"
fi

echo "" | tee -a "${log_file}"
echo "Log file: ${log_file}" | tee -a "${log_file}"
