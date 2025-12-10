#!/bin/bash
set -euo pipefail

## 07_render_ehh_report.sh
## Render EHH analysis R Markdown notebook to HTML
##
## Usage:
##   bash code/12_Sweep_Detection/ehh/07_render_ehh_report.sh
##
## Prerequisites:
##   - All previous steps completed (01-06)
##   - R with rmarkdown, tidyverse, patchwork packages installed
##
## Author: Jason Gallant Lab
## Date: 2024

root="$(git rev-parse --show-toplevel)"
source "${root}/pkings_gwas.env"

echo "=========================================="
echo "Rendering EHH Analysis Report"
echo "=========================================="
echo ""

# Input Rmd file
rmd_file="${root}/code/12_Sweep_Detection/ehh/ehh_analysis.Rmd"

# Output directory
output_dir="${root}/output_data/12_Sweep_Detection/ehh"
mkdir -p "${output_dir}"

# Check if Rmd file exists
if [[ ! -f "${rmd_file}" ]]; then
  echo "ERROR: R Markdown file not found: ${rmd_file}"
  exit 1
fi

echo "Input: ${rmd_file}"
echo "Output directory: ${output_dir}"
echo ""

# Render with R
echo "Rendering R Markdown..."
singularity exec --bind ${root}:/project_root ${rehh_image} \
  Rscript -e "rmarkdown::render(
  '/project_root/code/12_Sweep_Detection/ehh/ehh_analysis.Rmd',
  output_dir = '/project_root/output_data/12_Sweep_Detection/ehh',
  output_file = 'ehh_analysis.html'
)"

exit_code=$?

echo ""
if [[ ${exit_code} -eq 0 ]]; then
  echo "=========================================="
  echo "Report rendered successfully"
  echo "=========================================="
  echo ""
  echo "Output: ${output_dir}/ehh_analysis.html"
  echo ""
  echo "To view the report, open in a web browser:"
  echo "  open ${output_dir}/ehh_analysis.html"
else
  echo "=========================================="
  echo "ERROR: Report rendering failed"
  echo "=========================================="
  exit ${exit_code}
fi
