## region_stats_lib.R
## The search-region operator used for both the GWAS peaks and the pseudo-peaks.
##
## DESIGN.md section 5 requires that the peaks and the null be reduced by the *identical*
## operator: same region half-width, same window size, same step, same minimum number of scored
## sites, and the same aggregation over lineages. Everything here therefore takes those settings
## as arguments, and the callers pass one shared settings object. assert_same_operator() is
## called before any percentile is computed.
##
## A "statistic source" is a genome-wide table of (chrom, pos, value) for one lineage, plus:
##   direction  "high" if sweep-like means large, "low" if sweep-like means small
##   operator   "window_mean" (best W-bp sub-window of the per-site mean) or
##              "region_max"  (the extreme single grid value in the region)
##
## Loaded by 02_region_stats.R and reused by the calibration script.

suppressPackageStartupMessages({
  library(data.table)
})

## ----------------------------------------------------------------- operator settings
## One object, passed everywhere, so peaks and null cannot drift apart.
make_operator <- function(half_width = 100e3,
                          window     = 20e3,
                          step       =  5e3,
                          min_sites  =    20,
                          label      = "primary") {
  structure(list(half_width = half_width, window = window, step = step,
                 min_sites = min_sites, label = label),
            class = "region_operator")
}

format_operator <- function(op) {
  sprintf("%s: region +/-%g kb, window %g kb, step %g kb, min %d sites",
          op$label, op$half_width / 1e3, op$window / 1e3, op$step / 1e3, op$min_sites)
}

## Fail loudly if two result sets were produced by different operators.
assert_same_operator <- function(a, b, what = "peaks and null") {
  fields <- c("half_width", "window", "step", "min_sites")
  for (f in fields) {
    if (!isTRUE(all.equal(a[[f]], b[[f]]))) {
      stop(sprintf("operator mismatch between %s: %s is %s versus %s",
                   what, f, format(a[[f]]), format(b[[f]])), call. = FALSE)
    }
  }
  invisible(TRUE)
}

## ------------------------------------------------------------------- the operator itself
## For one region, reduce a (pos, value) table to a single number.
##
## window_mean: slide a window of op$window across the region in steps of op$step, take the mean
##   of the values whose positions fall inside it, require at least op$min_sites of them, and
##   return the largest (direction "high") or smallest (direction "low") window mean. Ties and
##   all-NA regions return NA.
## region_max: return the largest (or smallest) single value in the region. Used for grid
##   statistics such as SweepFinder2 CLR, where the grid is already a smoothed summary and a
##   further mean would double-smooth it.
##
## Returns a one-row list: value, the offset of the winning window centre from the region
## centre, and the number of scored sites in the region and in the winning window.
region_statistic <- function(pos, value, centre, op, direction = c("high", "low"),
                             operator = c("window_mean", "region_max")) {
  direction <- match.arg(direction)
  operator  <- match.arg(operator)
  lo <- centre - op$half_width
  hi <- centre + op$half_width

  keep <- which(pos >= lo & pos <= hi & !is.na(value))
  n_region <- length(keep)
  empty <- list(value = NA_real_, offset = NA_real_, n_region = n_region, n_window = NA_integer_)
  if (n_region < op$min_sites) return(empty)

  p <- pos[keep]; v <- value[keep]
  o <- order(p); p <- p[o]; v <- v[o]

  if (operator == "region_max") {
    i <- if (direction == "high") which.max(v) else which.min(v)
    return(list(value = v[i], offset = p[i] - centre,
                n_region = n_region, n_window = n_region))
  }

  ## sliding window means via cumulative sums: O(n + n_windows)
  starts <- seq(lo, hi - op$window, by = op$step)
  cs <- c(0, cumsum(v))
  a <- findInterval(starts - 1, p)              # sites strictly before the window
  b <- findInterval(starts + op$window - 1, p)  # sites up to the window end
  n <- b - a
  means <- ifelse(n >= op$min_sites, (cs[b + 1] - cs[a + 1]) / n, NA_real_)
  if (all(is.na(means))) return(empty)

  i <- if (direction == "high") which.max(means) else which.min(means)
  list(value = means[i], offset = starts[i] + op$window / 2 - centre,
       n_region = n_region, n_window = n[i])
}

## ---------------------------------------------------------- apply to many regions at once
## stat_dt: data.table with columns chrom, pos, value (one lineage, one statistic)
## regions: data.table with columns region_id, chrom, centre
apply_operator <- function(stat_dt, regions, op, direction, operator) {
  setkey(stat_dt, chrom, pos)
  out <- vector("list", nrow(regions))
  for (ch in unique(regions$chrom)) {
    s <- stat_dt[.(ch), nomatch = NULL]
    idx <- which(regions$chrom == ch)
    if (!nrow(s)) {
      for (i in idx) out[[i]] <- data.table(region_id = regions$region_id[i],
                                            value = NA_real_, offset = NA_real_,
                                            n_region = 0L, n_window = NA_integer_)
      next
    }
    for (i in idx) {
      r <- region_statistic(s$pos, s$value, regions$centre[i], op, direction, operator)
      out[[i]] <- data.table(region_id = regions$region_id[i], value = r$value,
                             offset = r$offset, n_region = r$n_region, n_window = r$n_window)
    }
  }
  rbindlist(out)
}

## --------------------------------------------------------------------- statistic loaders
## Each returns a data.table(chrom, pos, value) for one lineage.

load_rsb <- function(root, pair) {
  f <- file.path(root, "output_data/12_Sweep_Detection/ehh/merged",
                 sprintf("%s.genome_wide_rsb.txt.gz", pair))
  if (!file.exists(f)) stop("missing Rsb file: ", f)
  d <- fread(cmd = sprintf("zcat < %s", shQuote(f)))
  setnames(d, c("CHR", "POSITION", "RSB"), c("chrom", "pos", "value"), skip_absent = TRUE)
  d[!is.na(value), .(chrom, pos, value)]
}

load_clr <- function(root, lineage, tool = "sweepfinder2", suffix = "genome_wide_clr.txt") {
  f <- file.path(root, "output_data/12_Sweep_Detection", tool, "merged",
                 sprintf("%s.%s", lineage, suffix))
  if (!file.exists(f)) stop("missing CLR file: ", f)
  d <- fread(f)
  setnames(d, c("chromosome", "position", "CLR"), c("chrom", "pos", "value"), skip_absent = TRUE)
  d[!is.na(value), .(chrom, pos, value)]
}

## pi ratio from the lineage-group pixy run: BP pi / TP pi per 5-kb window, at window midpoints.
## Windows where either population has no callable sites are dropped rather than set to zero.
load_pi_ratio <- function(root, bp, tp, pixy_dir = "output_data/12_Sweep_Detection/pixy_lineages") {
  files <- list.files(file.path(root, pixy_dir), pattern = "_pi\\.txt$", full.names = TRUE)
  if (!length(files)) stop("no pixy pi files in ", file.path(root, pixy_dir))
  d <- rbindlist(lapply(files, fread), use.names = TRUE)
  d <- d[pop %in% c(bp, tp) & no_sites > 0 & !is.na(avg_pi)]
  w <- dcast(d, chromosome + window_pos_1 + window_pos_2 ~ pop,
             value.var = c("avg_pi", "no_sites"))
  bp_pi <- sprintf("avg_pi_%s", bp); tp_pi <- sprintf("avg_pi_%s", tp)
  w <- w[!is.na(get(bp_pi)) & !is.na(get(tp_pi)) & get(tp_pi) > 0]
  w[, .(chrom = chromosome,
        pos   = as.integer((window_pos_1 + window_pos_2) / 2),
        value = get(bp_pi) / get(tp_pi))]
}

## Callable sites per 5-kb window for one lineage, used as a matching covariate.
load_callable <- function(root, lineage, pixy_dir = "output_data/12_Sweep_Detection/pixy_lineages") {
  files <- list.files(file.path(root, pixy_dir), pattern = "_pi\\.txt$", full.names = TRUE)
  d <- rbindlist(lapply(files, fread), use.names = TRUE)[pop == lineage]
  d[, .(chrom = chromosome,
        pos   = as.integer((window_pos_1 + window_pos_2) / 2),
        value = as.numeric(no_sites))]
}

## -------------------------------------------------------------------------- test units
## The peaks, the lineages whose BP allele arose there, and each lineage's TP reference.
## DESIGN.md section 2.
test_units <- function() {
  data.table(
    peak = c(1L, 2L, 2L, 3L, 3L, 4L, 5L, 6L),
    bp   = c("BP1", "BP2", "BP3", "BP2", "BP3", "BP1", "BP1", "BP2"),
    tp   = c("TP1", "TP2", "TP2", "TP2", "TP2", "TP1", "TP1", "TP2")
  )
}

pair_of <- function(bp) c(BP1 = "BP1_TP1", BP2 = "BP2_TP2", BP3 = "BP3_TP2")[[bp]]

load_peaks <- function(root) {
  p <- fread(file.path(root, "input_data/10_Introgression/peaks.bed"),
             col.names = c("chrom_num", "start", "end", "peak", "score", "index_snp"))
  p[, chrom := paste0("chr", chrom_num)][order(peak)]
}

## ------------------------------------------------------------------------ region construction
## Peaks plus a deterministic grid of pseudo-peak centres. Shared by 02_region_stats.R and
## 02b_region_callable.R so that both score exactly the same regions; a mismatch would make the
## matching covariates line up with the wrong regions.
PSEUDO_SPACING_DEFAULT <- 200e3
PEAK_BUFFER_DEFAULT    <- 250e3

build_regions <- function(root, spacing = PSEUDO_SPACING_DEFAULT, buffer = PEAK_BUFFER_DEFAULT) {
  peaks <- load_peaks(root)
  fai <- fread(file.path(root, "input_data/00_Reference_Genome",
                         "jordan-mic4273-mb-hirise-20xx3__06-13-2024__final_assembly.fasta.renamed.fa.fai"),
               select = 1:2, col.names = c("chrom", "length"))
  fai <- fai[chrom %in% paste0("chr", 1:25)]

  peak_regions <- peaks[, .(region_id = paste0("peak", peak), region_type = "peak",
                            peak = peak, chrom = chrom, centre = as.numeric(index_snp))]

  grid <- fai[, .(centre = seq(spacing, length - spacing, by = spacing)), by = chrom]
  grid[, blocked := FALSE]
  for (i in seq_len(nrow(peaks))) {
    grid[chrom == peaks$chrom[i] &
           centre > peaks$start[i] - buffer &
           centre < peaks$end[i]   + buffer, blocked := TRUE]
  }
  pseudo_regions <- grid[blocked == FALSE,
                         .(region_id = paste0("pseudo_", chrom, "_", centre),
                           region_type = "pseudo", peak = NA_integer_, chrom, centre)]
  list(regions = rbind(peak_regions, pseudo_regions),
       n_blocked = sum(grid$blocked), peaks = peaks)
}
