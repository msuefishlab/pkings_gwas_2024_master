#!/bin/bash
set -euo pipefail

## 02_calculate_background_sfs.sh
## Calculate genome-wide background SFS for each population
## Uses SweepFinder2 -f option to compute empirical frequency spectrum
##
## This creates the null model (expected SFS) for genome-wide CLR calculations.
## The background SFS is calculated from all chromosome SFS files per population.
##
## Usage:
##   bash code/12_Sweep_Detection/02_calculate_background_sfs.sh
##
## Prerequisites:
##   - 01_convert_vcf_to_sf2.sh must be run first
##   - SweepFinder2 container must be built
##
## Author: Jason Gallant Lab
## Date: 2024

root="$(git rev-parse --show-toplevel)"
source "${root}/pkings_gwas.env"

# Input and output directories
indir="${root}/output_data/12_Sweep_Detection/sweepfinder2/sfs_input"
outdir="${root}/output_data/12_Sweep_Detection/sweepfinder2/background"
mkdir -p "${outdir}"

# Check if input directory exists and has files
if [[ ! -d "${indir}" ]]; then
  echo "ERROR: Input directory not found: ${indir}"
  echo "Please run 01_convert_vcf_to_sf2.sh first"
  exit 1
fi

# Check if container exists
if [[ ! -f "${sweepfinder2_image}" ]]; then
  echo "ERROR: SweepFinder2 container not found: ${sweepfinder2_image}"
  echo "Please build the container first:"
  echo "  singularity build ${sweepfinder2_image} images/sweepfinder2.def"
  exit 1
fi

# Population sets
POPS=(BP1 TP1 BP2 TP2 BP3)

# Log file
log_file="${outdir}/background_sfs.log"
echo "Starting background SFS calculation: $(date)" > "${log_file}"

# Process each population
for pop in "${POPS[@]}"; do
  echo "" | tee -a "${log_file}"
  echo "========================================" | tee -a "${log_file}"
  echo "Calculating background SFS for ${pop}" | tee -a "${log_file}"
  echo "========================================" | tee -a "${log_file}"

  # Check if SFS files exist for this population
  sfs_files=(${indir}/${pop}.chr*.sf2.txt)
  if [[ ! -e "${sfs_files[0]}" ]]; then
    echo "ERROR: No SFS files found for ${pop}" | tee -a "${log_file}"
    echo "  Expected pattern: ${indir}/${pop}.chr*.sf2.txt" | tee -a "${log_file}"
    continue
  fi

  echo "Found ${#sfs_files[@]} chromosome SFS files for ${pop}" | tee -a "${log_file}"

  # Concatenate all chromosomes for this population
  all_chroms="${outdir}/${pop}.allchroms.sf2.txt"
  echo "Concatenating chromosome SFS files..." | tee -a "${log_file}"

  # Write header (from first file)
  head -n 1 "${sfs_files[0]}" > "${all_chroms}"

  # Append all chromosome data (skip headers)
  for chr_file in "${sfs_files[@]}"; do
    chr_name=$(basename "${chr_file}" .sf2.txt)
    echo "  Adding ${chr_name}" | tee -a "${log_file}"
    tail -n +2 "${chr_file}" >> "${all_chroms}"
  done

  # Count total sites
  total_sites=$(tail -n +2 "${all_chroms}" | wc -l)
  echo "Total sites in concatenated file: ${total_sites}" | tee -a "${log_file}"

  # Run SweepFinder2 to calculate spectrum
  # -f option: compute allele frequency spectrum from input
  spectrum_file="${outdir}/${pop}.spectrum"
  echo "Running SweepFinder2 -f to calculate empirical frequency spectrum..." | tee -a "${log_file}"

  singularity exec --bind ${root}:/project_root ${sweepfinder2_image} \
    SweepFinder2 -f "/project_root/output_data/12_Sweep_Detection/sweepfinder2/background/${pop}.allchroms.sf2.txt" \
                    "/project_root/output_data/12_Sweep_Detection/sweepfinder2/background/${pop}.spectrum" \
    2>&1 | tee -a "${log_file}"

  # Check if spectrum file was created
  if [[ -f "${spectrum_file}" ]]; then
    echo "Background spectrum created: ${spectrum_file}" | tee -a "${log_file}"
    echo "Spectrum file size: $(du -h ${spectrum_file} | cut -f1)" | tee -a "${log_file}"
    echo "Spectrum file lines: $(wc -l < ${spectrum_file})" | tee -a "${log_file}"
  else
    echo "ERROR: Background spectrum file not created!" | tee -a "${log_file}"
  fi
done

echo "" | tee -a "${log_file}"
echo "========================================" | tee -a "${log_file}"
echo "Background SFS calculation complete: $(date)" | tee -a "${log_file}"
echo "========================================" | tee -a "${log_file}"
echo "" | tee -a "${log_file}"
echo "Output directory: ${outdir}/" | tee -a "${log_file}"
echo "Log file: ${log_file}" | tee -a "${log_file}"
echo "" | tee -a "${log_file}"

# Summary
echo "Background SFS files created:" | tee -a "${log_file}"
for pop in "${POPS[@]}"; do
  if [[ -f "${outdir}/${pop}.spectrum" ]]; then
    echo "  ${pop}.spectrum: $(du -h ${outdir}/${pop}.spectrum | cut -f1)" | tee -a "${log_file}"
  else
    echo "  ${pop}.spectrum: NOT CREATED" | tee -a "${log_file}"
  fi
done

echo "" | tee -a "${log_file}"
echo "Next step: Run 03_submit_sweepfinder2.sh to submit genome-wide CLR scan jobs" | tee -a "${log_file}"
