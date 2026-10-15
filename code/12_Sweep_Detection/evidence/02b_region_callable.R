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
## Input:  output_data/12_Sweep_Detection/pixy_lineages/lineages_chr*_pi.txt
## Output: output_data/12_Sweep_Detection/evidence/region_callable.tsv

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
cat("\nWrote ", file.path(dir_ev, "region_callable.tsv"), "\n", sep = "")
