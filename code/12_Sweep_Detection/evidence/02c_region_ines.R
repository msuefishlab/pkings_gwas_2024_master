#!/usr/bin/env Rscript
## 02c_region_ines.R
## Mean inES per search region, per lineage: the recombination-proxy matching covariate of
## DESIGN.md section 5.
##
## Why this covariate exists. Rsb, XP-EHH and the diversity ratio all respond to the local
## recombination rate as well as to selection: a region of low recombination carries long
## haplotypes, reduced diversity and elevated differentiation without any selection having
## acted. There is no recombination map for this assembly, so the rate cannot be matched
## directly. Integrated haplotype length in a lineage that is NOT part of the contrast under
## test is the available stand-in: it responds to the same local recombination environment, but
## the allele being tested is not sweeping in it, so a sweep at the locus does not inflate it.
##
## This is a different control from the differentiation-matched null. Differentiation is itself
## raised by a sweep, so matching on it bounds the ascertainment effect but conflates selection
## with recombination. The proxy separates them: it asks whether a peak's haplotype evidence is
## unusual among regions whose haplotype structure is similar in the other lineages.
##
## The reduction is a plain mean over the whole region, not the best sub-window. A covariate
## should describe the region's background, and taking a maximum would make it chase the same
## extremes the test statistics are looking for. It is produced through make_operator() with
## window = 2 * half_width so that exactly one window spans the region, which keeps it on the
## same code path as every other reduction and under the same minimum-site rule.
##
## Usage (inside gwas_tools.sif, from the repo root):
##   Rscript code/12_Sweep_Detection/evidence/02c_region_ines.R
##
## Input:  output_data/12_Sweep_Detection/ehh/merged/{lineage}.genome_wide_ines.txt.gz
## Output: output_data/12_Sweep_Detection/evidence/region_ines.tsv

suppressPackageStartupMessages({
  library(data.table)
  library(rprojroot)
})

## .git is a directory in a clone and a file in a linked worktree, and rprojroot's
## has_file() criterion matches only the latter. Accept either, so the script runs
## wherever the repository is checked out.
root <- find_root(has_file(".git") | has_dir(".git"))
source(file.path(root, "code/12_Sweep_Detection/evidence/region_stats_lib.R"))

HALF_WIDTH <- 100e3    # must equal the primary operator's half-width
MIN_SITES  <- 20L      # as for Rsb, the statistic this covariate is meant to condition

dir_ev <- file.path(root, "output_data/12_Sweep_Detection/evidence")
dir.create(dir_ev, recursive = TRUE, showWarnings = FALSE)

regions <- build_regions(root)$regions
cat("Regions:", nrow(regions), "(", regions[region_type == "peak", .N], "peaks )\n")

## one window spanning the whole region: a plain region mean, same code path as the statistics
op <- make_operator(half_width = HALF_WIDTH, window = 2 * HALF_WIDTH,
                    step = 5e3, min_sites = MIN_SITES, label = "ines_covariate")
cat("Operator: ", format_operator(op), "\n\n", sep = "")

ines_dir <- rsb_dir(root)   # same preference order as Rsb: new scans, else the published ones
cat("Reading inES from", ines_dir, "\n")

lineages <- c("BP1", "BP2", "BP3", "TP1", "TP2")
out <- list()
for (lg in lineages) {
  f <- file.path(ines_dir, sprintf("%s.genome_wide_ines.txt.gz", lg))
  if (!file.exists(f)) { cat("[skip]", lg, "- no file:", f, "\n"); next }
  d <- fread(cmd = sprintf("zcat < %s", shQuote(f)), select = c("CHR", "POSITION", "inES"))
  setnames(d, c("chrom", "pos", "value"))
  d <- d[!is.na(value)]
  r <- apply_operator(d, regions, op, "high", "window_mean")
  r[, lineage := lg]
  out[[lg]] <- r[, .(region_id, lineage, mean_ines = value, n_region)]
  cat(sprintf("  %-4s %8d scored sites | %4d of %d peaks, %5d of %d pseudo-peaks reduced\n",
              lg, nrow(d),
              r[region_id %like% "^peak" & !is.na(value), .N], regions[region_type == "peak", .N],
              r[region_id %like% "^pseudo" & !is.na(value), .N],
              regions[region_type == "pseudo", .N]))
  rm(d); invisible(gc())
}
if (!length(out)) stop("no inES files could be read")

ines <- rbindlist(out)

## Rank-normalize within lineage before anything averages across lineages. Genome-wide median
## inES differs 2-4 fold between these lineages (BP1 8043, BP2 7871, BP3 3231, TP1 3747,
## TP2 1814), so a raw mean across lineages would be dominated by whichever happens to have the
## longest haplotypes rather than by the region. Ranks are taken over every region at once;
## with six peaks against ~4,383 pseudo-peaks the peaks cannot appreciably shift the scale.
ines[, ines_rank := 100 * frank(mean_ines, na.last = "keep", ties.method = "average") /
       sum(!is.na(mean_ines)), by = lineage]

f_out <- file.path(dir_ev, "region_ines.tsv")
fwrite(ines, f_out, sep = "\t")
cat("\nWrote ", f_out, " (", nrow(ines), " rows)\n", sep = "")

cat("\nPeak rows (rank is the within-lineage percentile of mean inES across all regions):\n")
print(dcast(ines[region_id %like% "^peak"], region_id ~ lineage, value.var = "ines_rank"))
cat("\nmean inES, same rows:\n")
print(dcast(ines[region_id %like% "^peak"], region_id ~ lineage, value.var = "mean_ines"))
