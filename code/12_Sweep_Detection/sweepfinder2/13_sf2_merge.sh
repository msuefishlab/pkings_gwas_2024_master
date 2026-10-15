#!/bin/bash
set -euo pipefail

## 13_sf2_merge.sh
## Concatenate the per-chromosome CLR scans into one genome-wide table per lineage, in the
## column layout the null engine expects (chromosome, position, CLR, alpha).
##
## Usage (from the repo root):
##   bash code/12_Sweep_Detection/sweepfinder2/13_sf2_merge.sh [unfolded|folded]

root="$(git rev-parse --show-toplevel)"
source "${root}/pkings_gwas.env"

mode="${1:-unfolded}"
indir="${root}/output_data/12_Sweep_Detection/sweepfinder2/clr_${mode}"
outdir="${root}/output_data/12_Sweep_Detection/sweepfinder2/merged"
mkdir -p "${outdir}"

suffix="genome_wide_clr.txt"
[[ "${mode}" == "folded" ]] && suffix="genome_wide_clr_folded.txt"

for pop in BP1 TP1 BP2 TP2 BP3; do
  out="${outdir}/${pop}.${suffix}"
  printf 'chromosome\tposition\tCLR\talpha\n' > "${out}"
  n_chr=0
  for i in $(seq 1 25); do
    f="${indir}/${pop}.chr${i}.clr.txt"
    [[ -f "${f}" ]] || { echo "WARNING: missing ${f}"; continue; }
    # SweepFinder2 writes: location LR alpha
    awk -v c="chr${i}" 'NR > 1 && $1 != "location" {print c "\t" $1 "\t" $2 "\t" $3}' "${f}" >> "${out}"
    n_chr=$((n_chr + 1))
  done
  echo "${pop}: $(($(wc -l < "${out}") - 1)) grid points from ${n_chr} chromosomes -> ${out}"
done
