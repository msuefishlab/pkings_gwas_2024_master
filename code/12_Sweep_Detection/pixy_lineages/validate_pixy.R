#!/usr/bin/env Rscript
## validate_pixy.R
## Check the lineage-group pixy run against the published per-population run.
##
## BP3 and BIR contain exactly the same 16 fish, so their pi must agree window by window. Any
## disagreement is attributable to the software, not to the samples: the published run used
## pixy 1.2.10.beta2 inside the original container, and this one uses the same pixy version in
## a container rebuilt in 2026 whose pandas moved from 1.5.3 to 2.2.2.
##
## The other lineage groups have no published counterpart (BAM holds 67 fish including 25
## wobbles, while BP1 holds 42), so BP3 is the only available check.
##
## Usage (inside gwas_tools.sif, from the repo root):
##   Rscript code/12_Sweep_Detection/pixy_lineages/validate_pixy.R

suppressPackageStartupMessages({
  library(data.table)
  library(rprojroot)
})

root <- find_root(".git")
new_dir <- file.path(root, "output_data/12_Sweep_Detection/pixy_lineages")
published <- file.path(root, "output_data/05_PopGen",
                       "MIC4273_HAP2_NEWREF_MSU618_all_fish_all_sites_pi.txt")

files <- list.files(new_dir, pattern = "_pi\\.txt$", full.names = TRUE)
if (!length(files)) stop("no lineage pixy output in ", new_dir)
cat("Lineage pixy files:", length(files), "\n")

new <- rbindlist(lapply(files, fread), use.names = TRUE)[pop == "BP3"]
cat("BP3 windows (new run):", nrow(new), "\n")
if (!file.exists(published)) stop("published pi table not found: ", published)

old <- fread(published)[pop == "BIR"]
cat("BIR windows (published):", nrow(old), "\n")

m <- merge(new[, .(chromosome, window_pos_1, pi_new = avg_pi, sites_new = no_sites)],
           old[, .(chromosome, window_pos_1, pi_old = avg_pi, sites_old = no_sites)],
           by = c("chromosome", "window_pos_1"))
cat("Windows in common:", nrow(m), "\n\n")
if (!nrow(m)) stop("no overlapping windows; check chromosome naming")

both <- m[!is.na(pi_new) & !is.na(pi_old)]
d <- both[, pi_new - pi_old]
cat("pi, BP3 (new) versus BIR (published), over", nrow(both), "windows with both defined:\n")
cat(sprintf("  identical windows: %d (%.4f%%)\n", sum(d == 0), 100 * mean(d == 0)))
cat(sprintf("  max |difference|:  %.3e\n", max(abs(d))))
cat(sprintf("  mean difference:   %+.3e\n", mean(d)))
cat(sprintf("  correlation:       %.10f\n", cor(both$pi_new, both$pi_old)))
cat(sprintf("  callable sites identical in %.4f%% of windows\n",
            100 * mean(both$sites_new == both$sites_old)))

tol <- 1e-12
if (max(abs(d)) < tol) {
  cat("\nPASS: the rebuilt container reproduces the published pi exactly.\n")
} else if (max(abs(d)) < 1e-6) {
  cat("\nPASS (to floating-point tolerance): max difference", max(abs(d)), "\n")
} else {
  cat("\nFAIL: pi differs beyond floating-point tolerance. Worst windows:\n")
  print(head(both[order(-abs(pi_new - pi_old))], 10))
}
