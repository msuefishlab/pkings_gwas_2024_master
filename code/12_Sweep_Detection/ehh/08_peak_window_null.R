#!/usr/bin/env Rscript
## 08_peak_window_null.R
## Rank each GWAS peak's Rsb / inES against equal-size genomic windows
##
## 06_extract_peak_ehh.R compares the MAXIMUM |Rsb| (or inES) over all SNPs in a peak
## window (266-3,000 SNPs) with the distribution of SINGLE-SNP values. A window maximum is
## almost always above the single-SNP 99th percentile, so those percentiles are not a test.
## This script ranks each peak window against windows of the same size elsewhere in the
## genome, using the same window statistic.
##
## Method (all parameters are set below):
##   1. Peak windows: the extended peaks, input_data/10_Introgression/peaks.bed (decision D-002).
##   2. Primary null: the genome is tiled into non-overlapping windows of the peak's length in
##      bp, starting at position 1 of each chromosome. Windows overlapping any extended peak
##      +/- BUFFER_BP are dropped, as are windows with fewer than MIN_SNPS scored SNPs.
##   3. Statistics per window: n SNPs, mean, max, min and max|.| of Rsb; max inES of the BP and
##      of the TP population. Percentile = 100 x fraction of null windows whose statistic is
##      <= the peak's (for min Rsb, the fraction whose min is >= the peak's, so a high value
##      means unusually negative Rsb, i.e. longer TP haplotypes).
##   4. Sensitivity null: non-overlapping blocks of the same NUMBER of consecutive scored SNPs
##      as the peak window, outside the same exclusions (columns *_snpblock).
##   5. Background exceedance: fraction of the null windows whose max|Rsb| exceeds the
##      genome-wide single-SNP 99th percentile of |Rsb| (two-sided), and whose max Rsb exceeds
##      the single-SNP 99th percentile of signed Rsb (one-sided, the threshold drawn in Fig. S7).
## No random sampling is used; every quantity is a deterministic tiling of the genome.
##
## Inputs:
##   output_data/12_Sweep_Detection/ehh/merged/{BP1_TP1,BP2_TP2,BP3_TP2}.genome_wide_rsb.txt.gz
##   output_data/12_Sweep_Detection/ehh/scans/{pair}.{BP,TP}.chr*.ines.txt.gz
##     (per-pair inES, so TP2 values match the SNP set of each comparison; the merged
##      TP2.genome_wide_ines file mixes the BP2_TP2 and BP3_TP2 scans)
##
## Outputs: output_data/12_Sweep_Detection/ehh/window_null/
##   peak_window_null.tsv           one row per peak x comparison (18 rows)
##   genome_wide_rsb_thresholds.tsv single-SNP Rsb percentiles and 300-kb exceedance per comparison
##   README.md                      column definitions
##
## Usage:
##   Rscript code/12_Sweep_Detection/ehh/08_peak_window_null.R
##
## Prerequisites:
##   - 05_merge_ehh_results.sh completed (merged Rsb files and per-pair scans exist)
##
## Author: Jason Gallant Lab
## Date: 2026

suppressPackageStartupMessages({
  library(data.table)
  library(rprojroot)
})

root <- find_root(".git/index")

cat("========================================\n")
cat("Peak windows vs equal-size genome windows\n")
cat("========================================\n\n")

# ----------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------
BUFFER_BP <- 250000   # null windows may not overlap any extended peak +/- this many bp
MIN_SNPS  <- 20       # null windows need at least this many scored SNPs
BG_WIN_BP <- 300000   # fixed window size for the genome-wide exceedance summary

comparisons <- c("BP1_TP1", "BP2_TP2", "BP3_TP2")
merged_dir  <- file.path(root, "output_data/12_Sweep_Detection/ehh/merged")
scan_dir    <- file.path(root, "output_data/12_Sweep_Detection/ehh/scans")
peaks_file  <- file.path(root, "input_data/10_Introgression/peaks.bed")
outdir      <- file.path(root, "output_data/12_Sweep_Detection/ehh/window_null")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

genes <- c(`1` = "c6orf132", `2` = "itih3", `3` = "cebpg/cebpa",
           `4` = "cntnap5b", `5` = "sptbn4b", `6` = "naa40")

# Which comparison tests which peak. Peaks 1, 4, 5 arose in BP1 (Bambomo); Peaks 2, 3 are
# shared by BP2 and BP3; the Peak 6 BP haplotype is BP2-specific (BP3 kept for completeness).
role_of <- function(peak, comp) {
  if (comp == "BP1_TP1" && peak %in% c(1, 4, 5)) return("own")
  if (comp == "BP2_TP2" && peak %in% c(2, 3, 6)) return("own")
  if (comp == "BP3_TP2" && peak %in% c(2, 3))    return("own")
  if (comp == "BP3_TP2" && peak == 6)            return("own_completeness")
  if (comp == "BP1_TP1" && peak == 2)            return("control")
  "off_target"
}

# ----------------------------------------------------------------------
# Load peaks
# ----------------------------------------------------------------------
if (!file.exists(peaks_file)) stop("ERROR: extended peaks BED not found: ", peaks_file)
peaks <- fread(peaks_file, col.names = c("chromosome", "start", "end", "peak", "max_log_p", "idx_bp"))
if (!grepl("^chr", peaks$chromosome[1])) peaks[, chromosome := paste0("chr", chromosome)]
setorder(peaks, peak)
cat("Loaded", nrow(peaks), "extended peaks from", basename(peaks_file), "\n")
print(peaks)
cat("\n")

excl <- peaks[, .(CHR = chromosome, s = start - BUFFER_BP, e = end + BUFFER_BP)]

# ----------------------------------------------------------------------
# Window helpers
# ----------------------------------------------------------------------

# Flag rows (CHR, POSITION) that fall in an excluded interval
flag_excluded <- function(chr, pos) {
  ex <- logical(length(pos))
  for (i in seq_len(nrow(excl))) ex <- ex | (chr == excl$CHR[i] & pos >= excl$s[i] & pos <= excl$e[i])
  ex
}

# Non-overlapping bp tiles of length L; returns per-window stats of `val`
# (d must have CHR, POSITION, val with NA already removed)
bp_windows <- function(d, L, rsb = TRUE) {
  d[, wid := (POSITION - 1L) %/% L]
  if (rsb) {
    out <- d[, .(n = .N, mean = mean(val), max = max(val), min = min(val), maxabs = max(abs(val))), by = .(CHR, wid)]
  } else {
    out <- d[, .(n = .N, max = max(val)), by = .(CHR, wid)]
  }
  d[, wid := NULL]
  out[, `:=`(s = wid * L + 1, e = (wid + 1) * L)]
  ex <- logical(nrow(out))
  for (i in seq_len(nrow(excl))) ex <- ex | (out$CHR == excl$CHR[i] & out$e >= excl$s[i] & out$s <= excl$e[i])
  out[!ex & n >= MIN_SNPS]
}

# Non-overlapping blocks of n consecutive scored SNPs, never spanning an excluded interval
snp_blocks <- function(d, n, rsb = TRUE) {
  keep <- !d$excluded
  seg <- rleid(d$CHR, d$excluded)[keep]
  x <- d$val[keep]
  b <- data.table(seg = seg, x = x)
  b[, blk := (seq_len(.N) - 1L) %/% n, by = seg]
  if (rsb) {
    out <- b[, .(k = .N, mean = mean(x), max = max(x), min = min(x), maxabs = max(abs(x))), by = .(seg, blk)]
  } else {
    out <- b[, .(k = .N, max = max(x)), by = .(seg, blk)]
  }
  out[k == n]
}

pct_le <- function(null, obs) 100 * mean(null <= obs)
pct_ge <- function(null, obs) 100 * mean(null >= obs)

read_ines_pair <- function(pair, side) {
  files <- list.files(scan_dir, pattern = paste0("^", pair, "\\.", side, "\\.chr[0-9]+\\.ines\\.txt\\.gz$"),
                      full.names = TRUE)
  if (length(files) == 0) stop("No inES scans found for ", pair, " ", side, " in ", scan_dir)
  d <- rbindlist(lapply(files, fread, select = c("CHR", "POSITION", "inES")))
  d <- d[!is.na(inES), .(CHR, POSITION, val = inES)]
  setorder(d, CHR, POSITION)
  d[, excluded := flag_excluded(CHR, POSITION)]
  d
}

# ----------------------------------------------------------------------
# Main loop over comparisons
# ----------------------------------------------------------------------
results    <- list()
thresholds <- list()

for (comp in comparisons) {
  cat("========================================\n")
  cat("Comparison:", comp, "\n")
  cat("========================================\n")

  rsb_file <- file.path(merged_dir, paste0(comp, ".genome_wide_rsb.txt.gz"))
  if (!file.exists(rsb_file)) stop("ERROR: Rsb file not found: ", rsb_file)
  rsb <- fread(rsb_file)
  rsb_col <- grep("^RSB", names(rsb), value = TRUE)[1]
  rsb <- rsb[!is.na(get(rsb_col)), .(CHR, POSITION, val = get(rsb_col))]
  setorder(rsb, CHR, POSITION)
  rsb[, excluded := flag_excluded(CHR, POSITION)]
  cat("  Rsb SNPs (non-NA):", nrow(rsb), "\n")

  # Genome-wide single-SNP thresholds (all SNPs, peaks included, as in 06_extract_peak_ehh.R)
  q_abs <- quantile(abs(rsb$val), c(0.95, 0.99), names = FALSE)
  q_sgn <- quantile(rsb$val, c(0.01, 0.05, 0.95, 0.99), names = FALSE)
  bg <- bp_windows(rsb, BG_WIN_BP)
  thresholds[[comp]] <- data.table(
    comparison = comp, n_snps = nrow(rsb),
    absRsb_p95 = q_abs[1], absRsb_p99 = q_abs[2],
    Rsb_p01 = q_sgn[1], Rsb_p05 = q_sgn[2], Rsb_p95 = q_sgn[3], Rsb_p99 = q_sgn[4],
    bg300kb_n_windows = nrow(bg),
    bg300kb_frac_maxabs_gt_absP99 = mean(bg$maxabs > q_abs[2]),
    bg300kb_frac_max_gt_Rsb_p99   = mean(bg$max > q_sgn[4]),
    bg300kb_frac_maxabs_gt_absP95 = mean(bg$maxabs > q_abs[1])
  )
  cat(sprintf("  |Rsb| p95 = %.3f, p99 = %.3f; signed Rsb p95 = %.3f, p99 = %.3f\n",
              q_abs[1], q_abs[2], q_sgn[3], q_sgn[4]))
  cat(sprintf("  300-kb windows (n = %d) with a SNP above |Rsb| p99: %.1f%%; above signed p99: %.1f%%\n",
              nrow(bg), 100 * mean(bg$maxabs > q_abs[2]), 100 * mean(bg$max > q_sgn[4])))

  bp_pop <- sub("_.*", "", comp)
  tp_pop <- sub(".*_", "", comp)
  cat("  Loading per-pair inES scans for", bp_pop, "and", tp_pop, "\n")
  ines <- list(BP = read_ines_pair(comp, "BP"), TP = read_ines_pair(comp, "TP"))

  for (i in seq_len(nrow(peaks))) {
    pk <- peaks[i]
    L <- pk$end - pk$start + 1
    win <- rsb[CHR == pk$chromosome & POSITION >= pk$start & POSITION <= pk$end]
    n <- nrow(win)
    obs <- list(mean = mean(win$val), max = max(win$val), min = min(win$val), maxabs = max(abs(win$val)))
    imax <- which.max(abs(win$val))

    nb <- bp_windows(rsb, L)
    ns <- snp_blocks(rsb, n)

    row <- data.table(
      peak = pk$peak, gene = genes[as.character(pk$peak)], comparison = comp,
      role = role_of(pk$peak, comp),
      chrom = pk$chromosome, start = pk$start, end = pk$end, length_bp = L,
      n_snps = n,
      mean_Rsb = obs$mean, max_Rsb = obs$max, min_Rsb = obs$min,
      maxabs_Rsb_signed = win$val[imax], pos_maxabs_Rsb = win$POSITION[imax],
      n_null_windows = nrow(nb),
      pct_mean_Rsb = pct_le(nb$mean, obs$mean),
      pct_max_Rsb = pct_le(nb$max, obs$max),
      pct_maxabs_Rsb = pct_le(nb$maxabs, obs$maxabs),
      pct_min_Rsb_TPside = pct_ge(nb$min, obs$min),
      peak_exceeds_absP99 = obs$maxabs > q_abs[2],
      bg_frac_maxabs_gt_absP99 = mean(nb$maxabs > q_abs[2]),
      bg_frac_max_gt_Rsb_p99 = mean(nb$max > q_sgn[4]),
      n_null_snpblocks = nrow(ns),
      pct_mean_Rsb_snpblock = pct_le(ns$mean, obs$mean),
      pct_max_Rsb_snpblock = pct_le(ns$max, obs$max),
      pct_maxabs_Rsb_snpblock = pct_le(ns$maxabs, obs$maxabs)
    )

    # inES: max in the window for BP and TP, against equal-size windows of the same population
    for (side in c("BP", "TP")) {
      d <- ines[[side]]
      w_i <- d[CHR == pk$chromosome & POSITION >= pk$start & POSITION <= pk$end]
      mx <- max(w_i$val)
      nb_i <- bp_windows(d, L, rsb = FALSE)
      ns_i <- snp_blocks(d, nrow(w_i), rsb = FALSE)
      pop <- if (side == "BP") bp_pop else tp_pop
      row[, (paste0(side, "_pop")) := pop]
      row[, (paste0(side, "_n_snps_inES")) := nrow(w_i)]
      row[, (paste0(side, "_max_inES")) := mx]
      row[, (paste0(side, "_pct_max_inES")) := pct_le(nb_i$max, mx)]
      row[, (paste0(side, "_n_null_windows_inES")) := nrow(nb_i)]
      row[, (paste0(side, "_pct_max_inES_snpblock")) := pct_le(ns_i$max, mx)]
    }

    cat(sprintf("  Peak %d (%s) n=%d: mean %.2f [%.1f], max %.2f [%.1f], max|Rsb| %.2f [%.1f]; inES %s %.1f, %s %.1f\n",
                pk$peak, row$role, n, obs$mean, row$pct_mean_Rsb, obs$max, row$pct_max_Rsb,
                row$maxabs_Rsb_signed, row$pct_maxabs_Rsb,
                bp_pop, row$BP_pct_max_inES, tp_pop, row$TP_pct_max_inES))
    results[[paste(comp, pk$peak)]] <- row
  }
  cat("\n")
}

# ----------------------------------------------------------------------
# Save
# ----------------------------------------------------------------------
res <- rbindlist(results)
role_order <- c("own", "own_completeness", "control", "off_target")
res[, role_rank := match(role, role_order)]
setorder(res, role_rank, peak, comparison)
res[, role_rank := NULL]
num_cols <- names(res)[vapply(res, is.double, logical(1))]
res[, (num_cols) := lapply(.SD, function(x) signif(x, 6)), .SDcols = num_cols]

thr <- rbindlist(thresholds)
thr_num <- names(thr)[vapply(thr, is.double, logical(1))]
thr[, (thr_num) := lapply(.SD, function(x) signif(x, 6)), .SDcols = thr_num]

fwrite(res, file.path(outdir, "peak_window_null.tsv"), sep = "\t")
fwrite(thr, file.path(outdir, "genome_wide_rsb_thresholds.tsv"), sep = "\t")

readme <- c(
  "# Peak windows vs equal-size genome windows (Rsb, inES)",
  "",
  "Made by `code/12_Sweep_Detection/ehh/08_peak_window_null.R`. Peaks are the extended windows in `input_data/10_Introgression/peaks.bed` (D-002); Rsb comes from `merged/*.genome_wide_rsb.txt.gz` and inES from the per-pair scans in `scans/`.",
  sprintf("- **Null windows:** the genome is tiled into non-overlapping windows of the peak's bp length, starting at position 1 of each chromosome. Windows overlapping any extended peak +/- %d kb, or holding fewer than %d scored SNPs, are dropped. No random sampling is used.", BUFFER_BP / 1000, MIN_SNPS),
  "- **pct_\\*:** 100 x the fraction of null windows whose statistic is <= the peak's. `pct_max_Rsb` uses signed max Rsb (BP side, positive = longer BP haplotypes). `pct_maxabs_Rsb` uses max |Rsb| (either side). `pct_min_Rsb_TPside` is the fraction whose min Rsb is >= the peak's (TP side). `BP_/TP_pct_max_inES` ranks max inES in the window for each population against that population's windows.",
  "- **\\*_snpblock:** the same percentiles against non-overlapping blocks of the same number of consecutive scored SNPs (sensitivity check).",
  "- **maxabs_Rsb_signed:** the Rsb value with the largest |Rsb| in the window, sign kept (negative = longer TP haplotypes).",
  "- **bg_frac_maxabs_gt_absP99 / bg_frac_max_gt_Rsb_p99:** fraction of that peak's null windows with at least one SNP above the genome-wide single-SNP 99th percentile of |Rsb|, or of signed Rsb. `peak_exceeds_absP99` says whether the peak window itself does.",
  "- **role:** `own` = the lineage in which the peak's BP allele arose (P1, P4, P5 = BP1; P2, P3 = BP2 and BP3; P6 = BP2); `own_completeness` = P6 in BP3; `control` = BP1 at P2; `off_target` = the rest.",
  "- **genome_wide_rsb_thresholds.tsv:** single-SNP percentiles of |Rsb| and signed Rsb per comparison, and the fraction of 300-kb background windows with a SNP above each 99th percentile.",
  "- **Audit check:** the 2026-10-06 claims audit (section 1) used this bp tiling with >= 20 SNPs, no buffer and the signed max. That reproduces its 'Pct vs equal-size windows' column exactly; with the 250-kb buffer `pct_max_Rsb` differs by <= 0.1. Its '44-62% of 300-kb windows' is the one-sided `bg300kb_frac_max_gt_Rsb_p99`; the two-sided |Rsb| version is higher."
)
writeLines(readme, file.path(outdir, "README.md"))

cat("Saved:\n")
cat("  ", file.path(outdir, "peak_window_null.tsv"), "\n")
cat("  ", file.path(outdir, "genome_wide_rsb_thresholds.tsv"), "\n")
cat("  ", file.path(outdir, "README.md"), "\n\n")

cat("Genome-wide Rsb thresholds:\n")
print(thr)
cat("\nOwn-lineage and control rows:\n")
print(res[role != "off_target", .(peak, gene, comparison, role, n_snps,
                                  mean_Rsb = round(mean_Rsb, 2), max_Rsb = round(max_Rsb, 2),
                                  maxabs = round(maxabs_Rsb_signed, 2),
                                  p_mean = round(pct_mean_Rsb, 1), p_max = round(pct_max_Rsb, 1),
                                  p_maxabs = round(pct_maxabs_Rsb, 1),
                                  p_BP_inES = round(BP_pct_max_inES, 1), p_TP_inES = round(TP_pct_max_inES, 1),
                                  bg_abs99 = round(bg_frac_maxabs_gt_absP99, 2))])

cat("\n========================================\n")
cat("Window null complete\n")
cat("========================================\n")
