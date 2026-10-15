#!/usr/bin/env Rscript
## 02b_region_callable.R
## Callable sites per search region, per lineage: the second matching covariate.
##
## A region where fewer sites could be genotyped yields fewer scored SNPs and noisier
## statistics, in the peaks and in the null alike. Matching on it (DESIGN.md section 5) stops
## the comparison from being between a well-covered peak and poorly covered background.
##
## The count comes from pixy's no_sites column, which is the number of sites with data in each
## 5-kb window for that population, computed from the all-sites VCF.
##
## Usage (inside gwas_tools.sif, from the repo root):
##   Rscript code/12_Sweep_Detection/evidence/02b_region_callable.R
##
## It also writes the BP-TP differentiation of each region (mean Weir-Cockerham Fst from the
## same pixy run). That is NOT used by the headline null: a sweep raises differentiation, so
## matching on it would draw the null from other sweeps and from low-recombination regions.
## It exists for the differentiation-matched sensitivity row of DESIGN.md section 7, which
## bounds how much of the signal the peaks' own ascertainment could explain.
##
## Input:  output_data/12_Sweep_Detection/pixy_lineages/lineages_chr*_{pi,fst}.txt
## Output: output_data/12_Sweep_Detection/evidence/region_callable.tsv
##         output_data/12_Sweep_Detection/evidence/region_differentiation.tsv

suppressPackageStartupMessages({
  library(data.table)
  library(rprojroot)
})

root <- find_root(".git")
source(file.path(root, "code/12_Sweep_Detection/evidence/region_stats_lib.R"))

HALF_WIDTH <- 100e3   # must equal the primary operator's half-width
dir_ev <- file.path(root, "output_data/12_Sweep_Detection/evidence")
dir.create(dir_ev, recursive = TRUE, showWarnings = FALSE)

pixy_dir <- file.path(root, "output_data/12_Sweep_Detection/pixy_lineages")
files <- list.files(pixy_dir, pattern = "_pi\\.txt$", full.names = TRUE)
if (!length(files)) stop("no pixy pi files in ", pixy_dir, " (run run_pixy_lineages.sb first)")
cat("Reading", length(files), "pixy pi files\n")

pi <- rbindlist(lapply(files, fread), use.names = TRUE)
setnames(pi, c("chromosome", "window_pos_1", "window_pos_2"), c("chrom", "start", "end"),
         skip_absent = TRUE)
pi <- pi[, .(chrom, start, end, pop, no_sites = as.numeric(no_sites))]
cat("Windows:", nrow(pi), "across", uniqueN(pi$chrom), "chromosomes,",
    uniqueN(pi$pop), "lineages\n")

regions <- build_regions(root)$regions
cat("Regions:", nrow(regions), "(", regions[region_type == "peak", .N], "peaks )\n")

## Sum no_sites over the 5-kb windows whose midpoint falls in the region.
pi[, mid := (start + end) / 2]
setkey(pi, pop, chrom, mid)

out <- list()
for (lg in sort(unique(pi$pop))) {
  for (ch in unique(regions$chrom)) {
    w <- pi[.(lg, ch), nomatch = NULL]
    if (!nrow(w)) next
    r <- regions[chrom == ch]
    lo <- r$centre - HALF_WIDTH; hi <- r$centre + HALF_WIDTH
    cs <- c(0, cumsum(w$no_sites))
    a <- findInterval(lo - 1e-9, w$mid); b <- findInterval(hi, w$mid)
    out[[length(out) + 1]] <- data.table(region_id = r$region_id, lineage = lg,
                                         callable = cs[b + 1] - cs[a + 1],
                                         n_windows = b - a)
  }
  cat("  ", lg, ": done\n", sep = "")
}

res <- rbindlist(out)
fwrite(res, file.path(dir_ev, "region_callable.tsv"), sep = "\t")

cat("\nCallable sites per region (median over regions):\n")
print(res[, .(median_callable = median(callable), median_windows = median(n_windows),
              n_regions = .N), by = lineage])
cat("\nPeaks:\n")
print(dcast(res[region_id %like% "^peak"], region_id ~ lineage, value.var = "callable"))
## ------------------------------------------------------- BP-TP differentiation per region
fst_files <- list.files(pixy_dir, pattern = "_fst\\.txt$", full.names = TRUE)
if (length(fst_files)) {
  fst <- rbindlist(lapply(fst_files, fread), use.names = TRUE)
  setnames(fst, c("chromosome", "window_pos_1", "window_pos_2"), c("chrom", "start", "end"),
           skip_absent = TRUE)
  fst[, mid := (start + end) / 2]
  pairs <- list(BP1 = "TP1", BP2 = "TP2", BP3 = "TP2")
  dout <- list()
  for (bp in names(pairs)) {
    tp <- pairs[[bp]]
    f <- fst[(pop1 == bp & pop2 == tp) | (pop1 == tp & pop2 == bp)]
    f <- f[!is.na(avg_wc_fst)]
    if (!nrow(f)) { cat("  no Fst windows for ", bp, " vs ", tp, "\n", sep = ""); next }
    setkey(f, chrom, mid)
    for (ch in unique(regions$chrom)) {
      w <- f[.(ch), nomatch = NULL]
      if (!nrow(w)) next
      r <- regions[chrom == ch]
      lo <- r$centre - HALF_WIDTH; hi <- r$centre + HALF_WIDTH
      cs <- c(0, cumsum(w$avg_wc_fst))
      a <- findInterval(lo - 1e-9, w$mid); b <- findInterval(hi, w$mid)
      n <- b - a
      dout[[length(dout) + 1]] <- data.table(
        region_id = r$region_id, lineage = bp,
        mean_fst = fifelse(n > 0, (cs[b + 1] - cs[a + 1]) / pmax(n, 1), NA_real_),
        n_fst_windows = n)
    }
    cat("  ", bp, " vs ", tp, ": done\n", sep = "")
  }
  if (length(dout)) {
    dres <- rbindlist(dout)
    fwrite(dres, file.path(dir_ev, "region_differentiation.tsv"), sep = "\t")
    cat("\nMean BP-TP Fst at the peaks:\n")
    print(dcast(dres[region_id %like% "^peak"], region_id ~ lineage, value.var = "mean_fst"))
  }
} else {
  cat("\nNo pixy Fst files found; the differentiation sensitivity row will be unavailable.\n")
}

cat("\nWrote ", file.path(dir_ev, "region_callable.tsv"), "\n", sep = "")
