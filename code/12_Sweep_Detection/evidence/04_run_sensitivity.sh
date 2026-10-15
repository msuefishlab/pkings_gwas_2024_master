#!/bin/bash
set -euo pipefail

## 04_run_sensitivity.sh
## Re-run the whole reduction under each alternative configuration of DESIGN.md section 7.
##
## Every row recomputes the peaks AND the null together, so each has its own internally
## consistent operator. The results are never combined by taking the best: that would give the
## peaks a best-of-several advantage the null does not have. They are reported side by side,
## and any peak whose call changes is named in the write-up.
##
## Rows:
##   (default)        index SNP +/-100 kb, best 20-kb window, matched null  [the headline]
##   _r50             +/-50 kb
##   _w40             40-kb window
##   _unmatched       every pseudo-peak, matched on nothing
##   _differentiation matched additionally on BP-TP allele-frequency difference
##   _recombination   matched additionally on the recombination proxy (inES of the lineages
##                    outside the contrast). Unlike _differentiation, which is raised by a
##                    sweep and so conflates selection with recombination, this row holds the
##                    recombination environment fixed and leaves selection free to show.
##   _rsbgw           full sample, Rsb standardized genome-wide instead of per chromosome
##   _n16             haplotype statistics from the n = 16 down-sampled scans
##
## _n16 must be read against _rsbgw, not against the headline. The only down-sampled Rsb is the
## one 10_genome_stats.R recomputes from the scan RDS files, and that is standardized across all
## 25 chromosomes at once, whereas the headline Rsb is standardized within each chromosome.
## Comparing _n16 with the headline would therefore change the sample size and the
## standardization together, and a flipped call could not be attributed to either. _rsbgw is the
## same genome-wide standardization at the full sample size, so _rsbgw versus _n16 isolates the
## sample size.
##
## Usage (from the repo root):
##   bash code/12_Sweep_Detection/evidence/04_run_sensitivity.sh

root="$(git rev-parse --show-toplevel)"
source "${root}/pkings_gwas.env"
cd "${root}"

ev=code/12_Sweep_Detection/evidence
sx() { singularity exec --bind /mnt/research,/mnt/ufs18,/mnt/gs21 "${gwas_tools_image}" "$@"; }

echo "=== primary configuration"
sx Rscript ${ev}/02_region_stats.R
sx Rscript ${ev}/03_evidence_matrix.R

echo "=== region +/-50 kb"
sx Rscript ${ev}/02_region_stats.R --half_width 50000 --suffix _r50
sx Rscript ${ev}/03_evidence_matrix.R --suffix _r50

echo "=== 40-kb window"
sx Rscript ${ev}/02_region_stats.R --window 40000 --suffix _w40
sx Rscript ${ev}/03_evidence_matrix.R --suffix _w40

echo "=== unmatched null"
sx Rscript ${ev}/03_evidence_matrix.R --null_mode unmatched --out_suffix _unmatched

echo "=== differentiation-matched null (sensitivity only; a sweep creates differentiation)"
sx Rscript ${ev}/03_evidence_matrix.R --null_mode differentiation --out_suffix _differentiation

echo "=== recombination-proxy-matched null"
if [[ -f output_data/12_Sweep_Detection/evidence/region_ines.tsv ]]; then
  sx Rscript ${ev}/03_evidence_matrix.R --null_mode recombination --out_suffix _recombination
else
  echo "    skipped: run ${ev}/02c_region_ines.R first"
fi

echo "=== n = 16 down-sampled haplotype statistics, with its control row"
if ls output_data/12_Sweep_Detection/ehh/rsb_check/*.n16.genome_wide.txt.gz >/dev/null 2>&1; then
  echo "--- control: full sample, Rsb standardized genome-wide"
  sx Rscript ${ev}/02_region_stats.R --rsb_source rsb_check --suffix _rsbgw
  sx Rscript ${ev}/03_evidence_matrix.R --suffix _rsbgw
  echo "--- n = 16"
  sx Rscript ${ev}/02_region_stats.R --subsample 16 --suffix _n16
  sx Rscript ${ev}/03_evidence_matrix.R --suffix _n16
else
  echo "    skipped: no .n16 statistics yet. Run"
  echo "      sbatch --export=root=\$(pwd),SUBSAMPLE=16 code/12_Sweep_Detection/ehh/14_ehh_scan_array.sb"
  echo "    then Rscript code/12_Sweep_Detection/ehh/10_genome_stats.R --subsample 16"
fi

echo
echo "=== calls across configurations"
sx Rscript -e '
suppressPackageStartupMessages(library(data.table))
d <- file.path("output_data/12_Sweep_Detection/evidence")
fs <- list.files(d, pattern = "^evidence_matrix.*\\.tsv$", full.names = TRUE)
all <- rbindlist(lapply(fs, function(f) {
  x <- fread(f)
  if (!"supported" %in% names(x)) return(NULL)
  x[primary == TRUE, .(config = sub("^evidence_matrix_?", "", tools::file_path_sans_ext(basename(f))),
                       peak, supported = supported[1], strong = strong[1],
                       empirical_p = empirical_p[1], tp_label = tp_label[1]), by = peak][, -1]
}), fill = TRUE)
if (nrow(all)) {
  all[config == "", config := "primary"]
  print(dcast(unique(all), peak ~ config, value.var = "supported"))
  cat("\nEmpirical p by configuration:\n")
  print(dcast(unique(all), peak ~ config, value.var = "empirical_p"))
}'
