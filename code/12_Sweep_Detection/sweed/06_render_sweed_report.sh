#!/bin/bash
set -euo pipefail

## 06_render_sweed_report.sh
## Render SweeD analysis R Markdown notebook to HTML
##
## Usage:
##   bash code/12_Sweep_Detection/sweed/06_render_sweed_report.sh
##
## Prerequisites:
##   - 05_extract_peak_clr.R completed
##   - sweed_analysis.Rmd exists
##
## Author: Jason Gallant Lab
## Date: 2024

root="$(git rev-parse --show-toplevel)"

rmd_file="${root}/code/12_Sweep_Detection/sweed/sweed_analysis.Rmd"
output_dir="${root}/output_data/12_Sweep_Detection/sweed"

# Check if Rmd file exists
if [[ ! -f "${rmd_file}" ]]; then
  echo "ERROR: R Markdown file not found: ${rmd_file}"
  exit 1
fi

# Create output directory
mkdir -p "${output_dir}"

echo "========================================="
echo "Rendering SweeD Analysis Report"
echo "========================================="
echo "Input:  ${rmd_file}"
echo "Output: ${output_dir}/sweed_analysis.html"
echo ""

# Render R Markdown to HTML
Rscript -e "rmarkdown::render(
  '${rmd_file}',
  output_dir = '${output_dir}',
  output_file = 'sweed_analysis.html'
)"

exit_code=$?

echo ""
echo "========================================="
if [[ ${exit_code} -eq 0 ]]; then
  echo "SUCCESS: Report rendered successfully"
  echo ""
  echo "Output: ${output_dir}/sweed_analysis.html"
  echo "File size: $(du -h ${output_dir}/sweed_analysis.html | cut -f1)"
  echo ""
  echo "Open report with:"
  echo "  open ${output_dir}/sweed_analysis.html"
else
  echo "ERROR: Rendering failed with exit code ${exit_code}"
  echo "Check R Markdown file for errors"
fi
echo "========================================="

exit ${exit_code}
