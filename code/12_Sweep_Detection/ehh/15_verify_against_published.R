#!/usr/bin/env Rscript
## 15_verify_against_published.R
## Check that the re-run rehh scans reproduce the published ones.
##
## Why this exists. On 2026-10-15 the first run of the scan array produced inES values that
## correlated 0.963 with the published scan and matched none of them exactly, although the
## input VCFs were byte-identical in content, the marker and haplotype counts agreed, and
## rehh was the same version (3.2.2, unchanged on CRAN between 2021 and 2026-01-30). Running
## the same script by hand afterwards reproduced the published values exactly. The scan array
## had been submitted with a dependency on the per-chromosome split array but was observed
## running alongside it, so the early tasks read VCFs that were still being written.
##
## A corrupted read of that kind is silent: the scan completes, the counts look right, and only
## a comparison against a known-good result exposes it. So the pipeline now makes that
## comparison and fails loudly.
##
## Usage (inside gwas_tools.sif, from the repo root):
##   Rscript code/12_Sweep_Detection/ehh/15_verify_against_published.R [chr6 chr13 ...]
##
## Exit status 0 if every comparison matches exactly, 1 otherwise.

suppressPackageStartupMessages({
  library(data.table)
  library(rprojroot)
})

root <- find_root(".git")
args <- commandArgs(trailingOnly = TRUE)
chroms <- if (length(args)) args else c("chr6", "chr13", "chr24")

d_new <- file.path(root, "output_data/12_Sweep_Detection/ehh/merged")
d_pub <- file.path(root, "output_data/12_Sweep_Detection/ehh/merged_published")
if (!dir.exists(d_pub)) stop("no published scans to compare against: ", d_pub)

read_one <- function(dir, file) {
  f <- file.path(dir, file)
  if (!file.exists(f)) return(NULL)
  x <- fread(cmd = sprintf("zcat < %s", shQuote(f)))
  setnames(x, 1:3, c("chrom", "pos", "value"))
  x[, .(chrom, pos, value)]
}

compare <- function(file, label) {
  a <- read_one(d_new, file); b <- read_one(d_pub, file)
  if (is.null(a) || is.null(b)) {
    cat(sprintf("  %-34s SKIP (missing on one side)\n", label)); return(NA)
  }
  m <- merge(a[chrom %in% chroms], b[chrom %in% chroms], by = c("chrom", "pos"),
             suffixes = c("_new", "_pub"))
  k <- m[!is.na(value_new) & !is.na(value_pub)]
  if (!nrow(k)) { cat(sprintf("  %-34s SKIP (no overlap)\n", label)); return(NA) }
  ident <- mean(k$value_new == k$value_pub)
  mx <- max(abs(k$value_new - k$value_pub))
  cat(sprintf("  %-34s n=%7d identical=%6.2f%% max|diff|=%.3g cor=%.6f\n",
              label, nrow(k), 100 * ident, mx, cor(k$value_new, k$value_pub)))
  isTRUE(all.equal(ident, 1)) && mx == 0
}

cat("Comparing the re-run scans against the published ones on: ",
    paste(chroms, collapse = ", "), "\n\n", sep = "")

results <- c(
  compare("BP1.genome_wide_ines.txt.gz", "inES BP1"),
  compare("BP2.genome_wide_ines.txt.gz", "inES BP2"),
  compare("BP3.genome_wide_ines.txt.gz", "inES BP3"),
  compare("BP1_TP1.genome_wide_rsb.txt.gz", "Rsb BP1_TP1"),
  compare("BP2_TP2.genome_wide_rsb.txt.gz", "Rsb BP2_TP2"),
  compare("BP3_TP2.genome_wide_rsb.txt.gz", "Rsb BP3_TP2")
)

tested <- results[!is.na(results)]
cat("\n")
if (!length(tested)) {
  cat("NOTHING COMPARED: no overlapping files were found.\n"); quit(status = 1)
}
if (all(tested)) {
  cat("PASS: the re-run reproduces the published scans exactly on these chromosomes.\n")
  quit(status = 0)
}
cat("FAIL: the re-run does not reproduce the published scans.\n",
    "      Do not use these scans. Check that the per-chromosome VCFs were complete and\n",
    "      unchanged before the scan array started, then re-run the scans.\n", sep = "")
quit(status = 1)
