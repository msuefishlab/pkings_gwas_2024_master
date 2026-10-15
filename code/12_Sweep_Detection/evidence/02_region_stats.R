#!/usr/bin/env Rscript
## 02_region_stats.R
## Reduce every genome-wide statistic to one number per search region, for the six GWAS peaks
## and for several thousand pseudo-peaks, using the identical operator for both (DESIGN.md
## section 5).
##
## Pseudo-peaks are a deterministic grid of centres every 200 kb across chr1-chr25 (so the
## +/-100 kb regions abut but do not overlap), minus every GWAS peak +/- 250 kb. A grid rather
## than a random sample keeps the null reproducible without a seed. Regions that fall in
## assembly gaps or unmappable sequence drop out through the operator's minimum-site rule;
## there is no gap annotation for this assembly, so that is the only filter.
##
## Output is deliberately long and unaggregated: one row per region x lineage x statistic.
## Matching, aggregation over lineages, percentiles and the decision rule all happen in
## 03_evidence_matrix.R, so this script can be re-run as further statistics arrive.
##
## Usage (inside gwas_tools.sif, from the repo root):
##   Rscript code/12_Sweep_Detection/evidence/02_region_stats.R [--stats rsb,clr,pi_ratio]
##
## Output: output_data/12_Sweep_Detection/evidence/region_stats.tsv.gz
##         output_data/12_Sweep_Detection/evidence/region_stats_operators.tsv

suppressPackageStartupMessages({
  library(data.table)
  library(rprojroot)
})

root <- find_root(".git")
source(file.path(root, "code/12_Sweep_Detection/evidence/region_stats_lib.R"))

## --stats       comma-separated statistics (default: all that are available)
## --half_width  search-region half-width in bp (default 100000; sensitivity: 50000)
## --window      sub-window size in bp       (default  20000; sensitivity: 40000)
## --suffix      appended to the output file names, so sensitivity runs do not overwrite
##               the primary configuration
args <- commandArgs(trailingOnly = TRUE)
arg_of <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}
want       <- strsplit(arg_of("--stats", "rsb,clr,pi_ratio,xpehh,ihs,nsl,clr_folded"), ",")[[1]]
HALF_WIDTH <- as.numeric(arg_of("--half_width", 100e3))
WINDOW     <- as.numeric(arg_of("--window", 20e3))
SUFFIX     <- arg_of("--suffix", "")

out_dir <- file.path(root, "output_data/12_Sweep_Detection/evidence")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

PEAK_BUFFER <- 250e3   # exclusion around every GWAS peak

## ------------------------------------------------------------------ statistic registry
## lineages: "bp" for contrast statistics (BP versus its own TP), "all" for population
##   statistics, which are also computed on TP1 and TP2 so that the lineage-specificity label of
##   DESIGN.md section 6 has something to compare against.
## direction: is a sweep expected to make this large ("high") or small ("low")?
## operator:  how the region is reduced (see region_stats_lib.R)
## min_sites: scored sites required in a window; pixy statistics are 5-kb bins, not sites
stat_registry <- list(
  rsb        = list(lineages = "bp", direction = "high", operator = "window_mean", min_sites = 20L, primary = TRUE,
                    load = function(bp, tp) load_rsb(root, pair_of(bp))),
  clr        = list(lineages = "all", direction = "high", operator = "region_max",  min_sites =  3L, primary = TRUE,
                    load = function(bp, tp) load_clr(root, bp)),
  pi_ratio   = list(lineages = "bp", direction = "low",  operator = "window_mean", min_sites =  4L, primary = TRUE,
                    load = function(bp, tp) load_pi_ratio(root, bp, tp)),
  xpehh      = list(lineages = "bp", direction = "high", operator = "window_mean", min_sites = 20L, primary = FALSE,
                    load = function(bp, tp) load_derived(root, "xpehh", bp)),
  ihs        = list(lineages = "all", direction = "high", operator = "window_mean", min_sites = 20L, primary = FALSE,
                    load = function(bp, tp) load_derived(root, "ihs", bp)),
  nsl        = list(lineages = "all", direction = "high", operator = "window_mean", min_sites = 20L, primary = FALSE,
                    load = function(bp, tp) load_derived(root, "nsl", bp)),
  clr_folded = list(lineages = "all", direction = "high", operator = "region_max",  min_sites =  3L, primary = FALSE,
                    load = function(bp, tp) load_clr(root, bp, suffix = "genome_wide_clr_folded.txt"))
)

## Statistics written by 10_genome_stats.R / 12_nsl.py as {stat}/{lineage}.genome_wide.txt.gz
load_derived <- function(root, stat, lineage) {
  f <- file.path(root, "output_data/12_Sweep_Detection/ehh", stat,
                 sprintf("%s.genome_wide.txt.gz", lineage))
  if (!file.exists(f)) stop("missing ", stat, " file: ", f)
  d <- fread(cmd = sprintf("zcat < %s", shQuote(f)))
  setnames(d, names(d)[1:3], c("chrom", "pos", "value"))
  d[!is.na(value), .(chrom, pos, value)]
}

## ------------------------------------------------------------------------- the regions
peaks <- load_peaks(root)
units <- test_units()
lineages <- sort(unique(units$bp))
cat("Peaks:\n"); print(peaks[, .(peak, chrom, index_snp)])

## pseudo-peak spacing tracks the region width so that the regions abut without overlapping
PSEUDO_SPACING <- 2 * HALF_WIDTH
built <- build_regions(root, spacing = PSEUDO_SPACING, buffer = PEAK_BUFFER)
regions <- built$regions
peak_regions   <- regions[region_type == "peak"]
pseudo_regions <- regions[region_type == "pseudo"]
cat(sprintf("\nPseudo-peaks: %d candidates on %d chromosomes (%d blocked by peak +/-%g kb)\n",
            nrow(pseudo_regions), uniqueN(pseudo_regions$chrom), built$n_blocked, PEAK_BUFFER / 1e3))

## --------------------------------------------------------------------------- computation
op_rows <- list(); results <- list()

for (stat in want) {
  if (is.null(stat_registry[[stat]])) { cat("\n[skip] unknown statistic:", stat, "\n"); next }
  spec <- stat_registry[[stat]]
  op <- make_operator(half_width = HALF_WIDTH, window = WINDOW,
                      min_sites = spec$min_sites, label = stat)

  use <- if (identical(spec$lineages, "all")) c(lineages, sort(unique(units$tp))) else lineages
  for (bp in use) {
    tp <- units$tp[units$bp == bp][1]   # each BP lineage has one TP reference; NA for a TP
    dt <- tryCatch(spec$load(bp, tp), error = function(e) {
      cat("\n[skip]", stat, bp, "-", conditionMessage(e), "\n"); NULL })
    if (is.null(dt) || !nrow(dt)) next

    cat(sprintf("\n%-10s %-4s %s | %d scored sites\n", stat, bp, format_operator(op), nrow(dt)))
    r <- apply_operator(dt, regions, op, spec$direction, spec$operator)
    r[, `:=`(statistic = stat, lineage = bp, direction = spec$direction,
             operator = spec$operator, primary = spec$primary)]
    results[[length(results) + 1]] <- r
    op_rows[[length(op_rows) + 1]] <- data.table(
      statistic = stat, lineage = bp, half_width = op$half_width, window = op$window,
      step = op$step, min_sites = op$min_sites, direction = spec$direction,
      operator = spec$operator, n_sites_genome = nrow(dt))
    cat(sprintf("  peaks scored: %d of %d; pseudo-peaks scored: %d of %d\n",
                r[region_id %like% "^peak" & !is.na(value), .N], nrow(peak_regions),
                r[region_id %like% "^pseudo" & !is.na(value), .N], nrow(pseudo_regions)))
  }
}

if (!length(results)) stop("no statistics could be computed")

res <- merge(rbindlist(results), regions, by = "region_id")
setcolorder(res, c("region_id", "region_type", "peak", "chrom", "centre",
                   "statistic", "lineage", "value", "offset", "n_region", "n_window",
                   "direction", "operator", "primary"))
setorder(res, statistic, lineage, region_type, chrom, centre)

f_out <- file.path(out_dir, sprintf("region_stats%s.tsv.gz", SUFFIX))
fwrite(res, f_out, sep = "\t")
fwrite(rbindlist(op_rows), file.path(out_dir, sprintf("region_stats_operators%s.tsv", SUFFIX)), sep = "\t")

cat("\n", strrep("-", 70), "\nWrote ", f_out, " (", nrow(res), " rows)\n", sep = "")
cat("Statistics x lineages computed:\n")
print(dcast(res[!is.na(value), .N, by = .(statistic, lineage)], statistic ~ lineage,
            value.var = "N", fill = 0))
