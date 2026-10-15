#!/bin/bash
set -uo pipefail

## 00_run_all_downstream.sh
## Everything after the SLURM arrays: merge, derive the remaining statistics, reduce to
## regions, build the evidence matrix, and run the sensitivity rows.
##
## Each stage is skipped when its inputs are missing, and reports why, so this can be run
## before every array has finished and run again later without redoing completed work.
##
## Usage (from the repo root):
##   bash code/12_Sweep_Detection/evidence/00_run_all_downstream.sh [--wait]
##
##   --wait  block until no ehh_scan or sf2_scan job is left in the queue
##
## Run this on a dev node (ssh dev-amd20) or inside a SLURM job. The login gateway has no
## singularity, and every container stage would fail there one after another.

root="$(git rev-parse --show-toplevel)"
source "${root}/pkings_gwas.env"
cd "${root}"

ev=code/12_Sweep_Detection/evidence
D=output_data/12_Sweep_Detection
sx()  { singularity exec --bind /mnt/research,/mnt/ufs18,/mnt/gs21 "${gwas_tools_image}" "$@"; }
sxr() { singularity exec --bind /mnt/research,/mnt/ufs18,/mnt/gs21 "${rehh_image}" "$@"; }
step() { echo; echo "=== $* ==="; }

if ! command -v singularity >/dev/null 2>&1; then
  echo "ERROR: singularity is not available here." >&2
  echo "       The login gateway cannot run containers; use a dev node (ssh dev-amd20)" >&2
  echo "       or submit this as a SLURM job." >&2
  exit 1
fi

if [[ "${1:-}" == "--wait" ]]; then
  step "waiting for the ehh_scan and sf2_scan arrays"
  while [[ $(squeue -u "$USER" -h -n ehh_scan,sf2_scan -o '%i' | wc -l) -gt 0 ]]; do sleep 120; done
  echo "queue clear: $(date)"
fi

step "SweepFinder2: merge the per-chromosome scans"
n_clr=$(ls ${D}/sweepfinder2/clr_unfolded/*.clr.txt 2>/dev/null | wc -l)
if [[ "${n_clr}" -ge 125 ]]; then
  bash code/12_Sweep_Detection/sweepfinder2/13_sf2_merge.sh unfolded
else
  echo "skipped: ${n_clr} of 125 CLR files present"
fi
if [[ $(ls ${D}/sweepfinder2/clr_folded/*.clr.txt 2>/dev/null | wc -l) -ge 125 ]]; then
  bash code/12_Sweep_Detection/sweepfinder2/13_sf2_merge.sh folded
fi

step "rehh: merge inES and Rsb from the new scans"
## Down-sampled scans carry an .nN tag and are not merged: the merge produces the genome-wide
## tables the headline analysis reads, and those must come from the full sample.
n_rsb=$(ls ${D}/ehh/scans/*.rsb.txt.gz 2>/dev/null | grep -vc '\.n[0-9]\+\.' || true)
if [[ "${n_rsb}" -ge 75 ]]; then
  bash code/12_Sweep_Detection/ehh/05_merge_ehh_results.sh || echo "merge reported a problem; continuing"
else
  echo "skipped: ${n_rsb} of 75 full-sample Rsb files present"
fi

step "verify the re-run scans against the published ones"
if [[ -d ${D}/ehh/merged_published ]] && ls ${D}/ehh/merged/*genome_wide_rsb* >/dev/null 2>&1; then
  if ! sx Rscript code/12_Sweep_Detection/ehh/15_verify_against_published.R; then
    echo "STOPPING: the re-run scans do not reproduce the published ones." >&2
    exit 1
  fi
else
  echo "skipped: nothing to compare yet"
fi

step "rehh: real iHS, XP-EHH and the Rsb reproduction check"
if [[ $(ls ${sweep_scratch}/12_Sweep_Detection/ehh_rds/*.scan.rds 2>/dev/null | grep -vc '\.n16\.') -ge 150 ]]; then
  sxr Rscript code/12_Sweep_Detection/ehh/10_genome_stats.R
else
  echo "skipped: full-sample scan RDS files not all present"
fi
if [[ $(ls ${sweep_scratch}/12_Sweep_Detection/ehh_rds/*.n16.*.scan.rds 2>/dev/null | wc -l) -ge 150 ]]; then
  sxr Rscript code/12_Sweep_Detection/ehh/10_genome_stats.R --subsample 16
fi

step "allele-specific EHH for the figure"
if [[ -f ${D}/by_chrom/BP1.polarized.chr6.vcf.gz ]]; then
  sxr Rscript code/12_Sweep_Detection/ehh/11_allele_specific_ehh.R || echo "allele-specific EHH failed; continuing"
fi

step "index-SNP facts and the applicability map"
sx Rscript ${ev}/01_index_snp_facts.R || echo "index-SNP facts failed; continuing"

step "matching covariates"
sx Rscript ${ev}/02b_region_callable.R || echo "covariates failed; continuing"

step "region statistics and the evidence matrix"
sx Rscript ${ev}/02_region_stats.R || { echo "region statistics FAILED"; exit 1; }
sx Rscript ${ev}/03_evidence_matrix.R || { echo "evidence matrix FAILED"; exit 1; }

step "sensitivity rows"
bash ${ev}/04_run_sensitivity.sh || echo "some sensitivity rows failed"

step "done: $(date)"
ls -la ${D}/evidence/
