#!/usr/bin/env Rscript
## 03_evidence_matrix.R
## Turn the per-region statistics into the evidence matrix, and measure how often the same
## decision rule fires on genome regions that are not GWAS peaks (DESIGN.md sections 5 and 6).
##
## The calibration is the point of this script. A percentile only means something once we know
## how often a non-peak region of the same kind reaches it, so every number reported for the
## peaks is also computed for several thousand pseudo-peaks, through the identical operator.
##
## Matching (DESIGN.md section 5): the headline null keeps only pseudo-peaks in the same decile
## as the peak for SNP density and callable fraction. It is deliberately NOT matched on BP-TP
## differentiation: a sweep creates differentiation, so matching on it would build the null out
## of other sweeps. The differentiation-matched null is reported as a sensitivity row instead.
##
## Usage (inside gwas_tools.sif, from the repo root):
##   Rscript code/12_Sweep_Detection/evidence/03_evidence_matrix.R
##
## Input:  output_data/12_Sweep_Detection/evidence/region_stats.tsv.gz
## Output: evidence_matrix.tsv, calibration.tsv, statistic_correlations.tsv

suppressPackageStartupMessages({
  library(data.table)
  library(rprojroot)
})

## .git is a directory in a clone and a file in a linked worktree, and rprojroot's
## has_file() criterion matches only the latter. Accept either, so the script runs
## wherever the repository is checked out.
root <- find_root(has_file(".git") | has_dir(".git"))
source(file.path(root, "code/12_Sweep_Detection/evidence/region_stats_lib.R"))

## --suffix     read region_stats<suffix>.tsv.gz (the operator configuration)
## --out_suffix appended to the output names; defaults to --suffix. Give sensitivity runs
##              their own, or they overwrite the primary configuration's results.
## --null_mode  "matched" (default; SNP density and callable fraction),
##              "unmatched" (every pseudo-peak),
##              "differentiation" (adds BP-TP Fst to the matching; a sensitivity row only,
##              because a sweep creates differentiation), or
##              "recombination" (adds the recombination proxy of DESIGN.md section 5: mean
##              rank-normalized inES over the lineages not in the peak's contrast, from
##              02c_region_ines.R). The proxy responds to the local recombination environment
##              but not to a sweep on the allele under test, so unlike the differentiation
##              row it separates those two explanations instead of conflating them.
args <- commandArgs(trailingOnly = TRUE)
arg_of <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}
SUFFIX     <- arg_of("--suffix", "")
OUT_SUFFIX <- arg_of("--out_suffix", SUFFIX)
NULL_MODE  <- arg_of("--null_mode", "matched")

dir_ev  <- file.path(root, "output_data/12_Sweep_Detection/evidence")
stats   <- fread(file.path(dir_ev, sprintf("region_stats%s.tsv.gz", SUFFIX)))
units   <- test_units()
cat("Region statistics: region_stats", SUFFIX, ".tsv.gz | null mode: ", NULL_MODE,
    " | output suffix: '", OUT_SUFFIX, "'\n", sep = "")

PRIMARY      <- c("rsb", "clr", "pi_ratio")   # one per evidence class (H, S, D)
SUPPORT_PCT  <- 95                            # "supported": >= 2 of 3 primaries at this percentile
STRONG_PCT   <- 99                            # "strong": supported and >= 1 primary here
N_DECILES    <- 10

## -------------------------------------------------------------- percentile within a null set
## Fraction of null values at or below the observed one, in the direction that means sweep-like.
## For "low" statistics (the pi ratio) a small value is extreme, so the sense is reversed.
pct_of <- function(value, null_values, direction) {
  null_values <- null_values[!is.na(null_values)]
  if (is.na(value) || !length(null_values)) return(NA_real_)
  if (direction == "high") 100 * mean(null_values <= value)
  else                     100 * mean(null_values >= value)
}

## ------------------------------------------------------------- aggregate over eligible lineages
## A peak's statistic is the most sweep-like value over the lineages in which its BP allele
## arose; the null must take the same maximum over the same lineages, or the peak gets a
## best-of-several advantage the null does not have.
aggregate_over_lineages <- function(dt, lineage_set, direction) {
  x <- dt[lineage %in% lineage_set]
  if (!nrow(x)) return(data.table(region_id = character(), value = numeric()))
  if (direction == "high") x[, .(value = max(value, na.rm = TRUE)), by = region_id]
  else                     x[, .(value = min(value, na.rm = TRUE)), by = region_id]
}

## ------------------------------------------------------------------------ matching covariates
## SNP density comes from the operator's own site count, so it is on exactly the scale that
## affects the statistic. Callable fraction is the number of sites pixy could score.
covars <- stats[statistic == "rsb", .(n_snps = mean(n_region, na.rm = TRUE)),
                by = .(region_id, region_type, peak, chrom, centre, lineage)]
if (file.exists(file.path(dir_ev, "region_callable.tsv"))) {
  callable <- fread(file.path(dir_ev, "region_callable.tsv"))
  covars <- merge(covars, callable, by = c("region_id", "lineage"), all.x = TRUE)
} else {
  cat("NOTE: region_callable.tsv not found; matching on SNP density only.\n",
      "      Run 02b_region_callable.R once the lineage pixy output exists.\n", sep = "")
  covars[, callable := NA_real_]
}
if (file.exists(file.path(dir_ev, "region_differentiation.tsv"))) {
  diffs <- fread(file.path(dir_ev, "region_differentiation.tsv"))
  covars <- merge(covars, diffs[, .(region_id, lineage, mean_fst)],
                  by = c("region_id", "lineage"), all.x = TRUE)
} else {
  covars[, mean_fst := NA_real_]
}

## The recombination proxy is kept out of `covars` on purpose. Every other covariate is
## averaged over the lineages being tested; this one is averaged over the lineages that are
## NOT being tested, and that set differs from peak to peak, so it is assembled inside the
## per-peak loop instead.
if (file.exists(file.path(dir_ev, "region_ines.tsv"))) {
  ines <- fread(file.path(dir_ev, "region_ines.tsv"))
  ines_lineages <- sort(unique(ines$lineage))
  cat("Recombination proxy available for: ", paste(ines_lineages, collapse = ", "), "\n", sep = "")
} else {
  ines <- data.table(); ines_lineages <- character()
}

decile <- function(x) {
  if (all(is.na(x))) return(rep(NA_integer_, length(x)))
  cut(x, breaks = unique(quantile(x, probs = seq(0, 1, length.out = N_DECILES + 1), na.rm = TRUE)),
      include.lowest = TRUE, labels = FALSE)
}

## ============================================================================ per-peak testing
matrix_rows <- list(); null_scores <- list()

for (pk in sort(unique(units$peak))) {
  lineage_set <- units[peak == pk, bp]
  tp_set      <- unique(units[peak == pk, tp])
  cat("\n", strrep("=", 68), "\nPeak ", pk, " | lineage(s): ", paste(lineage_set, collapse = ", "),
      " | matched TP: ", paste(tp_set, collapse = ", "), "\n", sep = "")

  ## --- matching: deciles are defined on the pseudo-peaks, using the same lineage set
  cv <- covars[lineage %in% lineage_set,
               .(n_snps = mean(n_snps, na.rm = TRUE), callable = mean(callable, na.rm = TRUE),
                 mean_fst = mean(mean_fst, na.rm = TRUE)),
               by = .(region_id, region_type, peak)]
  ## Recombination proxy for this peak: the lineages outside its contrast. They share the
  ## region's recombination environment but are not sweeping for the allele under test, so the
  ## proxy tracks recombination without absorbing the signal. DESIGN.md section 5 forbids
  ## taking it from the matched TP, which may share the signal; setdiff() excludes both the
  ## eligible BP lineages and that TP.
  proxy_set <- setdiff(ines_lineages, c(lineage_set, tp_set))
  if (length(proxy_set) && nrow(ines)) {
    prox <- ines[lineage %in% proxy_set, .(recomb = mean(ines_rank, na.rm = TRUE)), by = region_id]
    cv <- merge(cv, prox, by = "region_id", all.x = TRUE)
  } else {
    cv[, recomb := NA_real_]
  }

  cvp <- cv[region_type == "pseudo"]
  cvp[, `:=`(d_snps = decile(n_snps), d_call = decile(callable), d_fst = decile(mean_fst),
             d_recomb = decile(recomb))]
  target <- cv[region_type == "peak" & peak == pk]
  if (!nrow(target)) { cat("  no covariates for this peak; skipped\n"); next }

  ## Caliper matching, widened until the null is large enough to resolve a 99th percentile.
  ## An exact single-decile match leaves only a few hundred regions, so the 99th percentile
  ## would rest on four or five of them. Starting at the peak's own decile, the caliper grows
  ## by one decile at a time until at least MIN_NULL regions qualify; the width actually used
  ## is recorded, because a wide caliper means the matching did less work.
  MIN_NULL <- 1000
  dec_of <- function(x, ref) {
    br <- quantile(ref, seq(0, 1, length.out = N_DECILES + 1), na.rm = TRUE)[-c(1, N_DECILES + 1)]
    findInterval(x, br) + 1L
  }
  d_snps_t <- dec_of(target$n_snps, cvp$n_snps)
  use_call <- !all(is.na(cvp$callable)) && !is.na(target$callable)
  d_call_t <- if (use_call) dec_of(target$callable, cvp$callable) else NA_integer_
  ## differentiation is added only for the sensitivity row, never for the headline
  use_fst  <- NULL_MODE == "differentiation" && !all(is.na(cvp$mean_fst)) && !is.na(target$mean_fst)
  d_fst_t  <- if (use_fst) dec_of(target$mean_fst, cvp$mean_fst) else NA_integer_
  ## the recombination proxy is likewise added only for its own sensitivity row
  use_rec  <- NULL_MODE == "recombination" && nrow(target) == 1L &&
              !all(is.na(cvp$recomb)) && !is.na(target$recomb)
  d_rec_t  <- if (use_rec) dec_of(target$recomb, cvp$recomb) else NA_integer_
  if (NULL_MODE == "recombination" && !use_rec) {
    stop("--null_mode recombination needs region_ines.tsv with a value for peak ", pk,
         "; run 02c_region_ines.R first", call. = FALSE)
  }

  if (NULL_MODE == "unmatched") {
    matched_ids <- cvp$region_id
    caliper <- NA_integer_
  } else {
  caliper <- 0L
  repeat {
    ok <- abs(cvp$d_snps - d_snps_t) <= caliper
    if (use_call) ok <- ok & abs(cvp$d_call - d_call_t) <= caliper
    if (use_fst)  ok <- ok & abs(cvp$d_fst  - d_fst_t)  <= caliper
    if (use_rec)  ok <- ok & abs(cvp$d_recomb - d_rec_t) <= caliper
    matched_ids <- cvp$region_id[which(ok)]
    if (length(matched_ids) >= MIN_NULL || caliper >= N_DECILES) break
    caliper <- caliper + 1L
  }
  }
  match_note <- if (NULL_MODE == "unmatched") "unmatched (all pseudo-peaks)" else
    sprintf("%s, caliper +/-%d decile(s)",
                        if (use_call) "SNP density and callable fraction" else
                          "SNP density only (callable fraction unavailable)",
                        caliper)
  if (use_fst) match_note <- paste(match_note, "+ BP-TP Fst")
  if (use_rec) match_note <- paste0(match_note, " + recombination proxy (inES of ",
                                    paste(proxy_set, collapse = "/"), ")")
  cat("  matched null: ", length(matched_ids), " of ", nrow(cvp),
      " pseudo-peaks (", match_note, ")\n", sep = "")

  ## --- percentiles for each statistic, peak and every pseudo-peak
  per_stat <- list()
  for (st in unique(stats$statistic)) {
    spec_dir <- stats[statistic == st, direction][1]
    agg <- aggregate_over_lineages(stats[statistic == st], lineage_set, spec_dir)
    if (!nrow(agg)) next
    agg <- merge(agg, unique(stats[, .(region_id, region_type)]), by = "region_id")
    peak_val <- agg[region_id == paste0("peak", pk), value]
    if (!length(peak_val)) next

    null_all     <- agg[region_type == "pseudo", value]
    null_matched <- agg[region_id %in% matched_ids, value]

    ## every pseudo-peak also needs its own percentile, so the rule can be applied to it
    null_pcts <- if (spec_dir == "high") {
      100 * frank(null_matched, na.last = "keep", ties.method = "max") / sum(!is.na(null_matched))
    } else {
      100 * frank(-null_matched, na.last = "keep", ties.method = "max") / sum(!is.na(null_matched))
    }

    per_stat[[st]] <- list(
      statistic = st, direction = spec_dir,
      peak_value = peak_val,
      pct_matched   = pct_of(peak_val, null_matched, spec_dir),
      pct_unmatched = pct_of(peak_val, null_all,     spec_dir),
      null_ids = agg[region_id %in% matched_ids, region_id], null_pcts = null_pcts)

    ## TP control: the same statistic in the matched TP population, where it is defined per
    ## population rather than per contrast (CLR, iHS, nSL)
    tp_vals <- stats[statistic == st & lineage %in% tp_set]
    if (nrow(tp_vals)) {
      tp_agg  <- aggregate_over_lineages(tp_vals, tp_set, spec_dir)
      tp_agg  <- merge(tp_agg, unique(stats[, .(region_id, region_type)]), by = "region_id")
      tp_peak <- tp_agg[region_id == paste0("peak", pk), value]
      per_stat[[st]]$tp_value <- if (length(tp_peak)) tp_peak else NA_real_
      per_stat[[st]]$tp_pct   <- pct_of(per_stat[[st]]$tp_value,
                                        tp_agg[region_id %in% matched_ids, value], spec_dir)
    } else {
      per_stat[[st]]$tp_value <- NA_real_; per_stat[[st]]$tp_pct <- NA_real_
    }
  }
  if (!length(per_stat)) { cat("  no statistics available for this peak\n"); next }

  for (st in names(per_stat)) {
    p <- per_stat[[st]]
    matrix_rows[[length(matrix_rows) + 1]] <- data.table(
      peak = pk, lineages = paste(lineage_set, collapse = "+"), statistic = st,
      primary = st %in% PRIMARY, direction = p$direction,
      peak_value = p$peak_value, pct_matched = p$pct_matched, pct_unmatched = p$pct_unmatched,
      tp_value = p$tp_value, tp_pct = p$tp_pct,
      n_matched_null = length(p$null_ids), match_note = match_note)
    cat(sprintf("  %-11s %-8s value %9.3f  pct(matched) %5.1f  pct(all) %5.1f  TP pct %s\n",
                st, if (st %in% PRIMARY) "primary" else "second.",
                p$peak_value, p$pct_matched, p$pct_unmatched,
                ifelse(is.na(p$tp_pct), "   -", sprintf("%5.1f", p$tp_pct))))
  }

  ## --- the decision rule, applied to the peak and to every matched pseudo-peak
  prim <- intersect(PRIMARY, names(per_stat))
  if (length(prim) < 2) {
    cat("  fewer than two primary statistics available; no call made yet\n")
    next
  }
  peak_pcts <- vapply(prim, function(s) per_stat[[s]]$pct_matched, numeric(1))
  composite <- sort(peak_pcts, decreasing = TRUE)[2]        # the value the "2 of 3" rule turns on
  supported <- sum(peak_pcts >= SUPPORT_PCT, na.rm = TRUE) >= 2
  strong    <- supported && any(peak_pcts >= STRONG_PCT, na.rm = TRUE)

  ids <- Reduce(intersect, lapply(prim, function(s) per_stat[[s]]$null_ids))
  nm  <- do.call(cbind, lapply(prim, function(s) {
    v <- per_stat[[s]]$null_pcts; names(v) <- per_stat[[s]]$null_ids; v[ids]
  }))
  colnames(nm) <- prim
  null_comp <- apply(nm, 1, function(r) if (all(is.na(r))) NA_real_ else sort(r, decreasing = TRUE)[2])
  null_sup  <- apply(nm, 1, function(r) sum(r >= SUPPORT_PCT, na.rm = TRUE) >= 2)
  null_str  <- null_sup & apply(nm, 1, function(r) any(r >= STRONG_PCT, na.rm = TRUE))
  emp_p     <- mean(null_comp >= composite, na.rm = TRUE)

  tp_high <- vapply(prim, function(s) isTRUE(per_stat[[s]]$tp_pct >= SUPPORT_PCT), logical(1))
  label   <- if (!supported) "unsupported" else if (any(tp_high)) "shared with TP" else "BP-specific"

  cat(sprintf("\n  primaries used: %s\n  composite (2nd highest percentile): %.1f\n",
              paste(prim, collapse = ", "), composite))
  cat(sprintf("  call: %s | %s | empirical p = %.4f (%d matched pseudo-peaks)\n",
              if (strong) "STRONG" else if (supported) "supported" else "not supported",
              label, emp_p, length(ids)))
  cat(sprintf("  this rule fires on %.1f%% of matched pseudo-peaks (supported), %.1f%% (strong)\n",
              100 * mean(null_sup, na.rm = TRUE), 100 * mean(null_str, na.rm = TRUE)))

  for (i in seq_along(matrix_rows)) if (matrix_rows[[i]]$peak[1] == pk) {
    matrix_rows[[i]][, `:=`(composite = composite, supported = supported, strong = strong,
                            empirical_p = emp_p, tp_label = label)]
  }
  null_scores[[as.character(pk)]] <- data.table(
    peak = pk, n_null = length(ids),
    fpr_supported = mean(null_sup, na.rm = TRUE), fpr_strong = mean(null_str, na.rm = TRUE),
    null_composite_median = median(null_comp, na.rm = TRUE))
}

## ================================================================================= outputs
if (length(matrix_rows)) {
  em <- rbindlist(matrix_rows, fill = TRUE)
  fwrite(em, file.path(dir_ev, sprintf("evidence_matrix%s.tsv", OUT_SUFFIX)), sep = "\t")
  cat("\n", strrep("=", 68), "\nEvidence matrix (primary statistics):\n", sep = "")
  if ("supported" %in% names(em)) {
    print(dcast(em[primary == TRUE], peak + lineages + supported + strong + empirical_p + tp_label ~
                  statistic, value.var = "pct_matched"))
  } else {
    ## no peak had two primary statistics available yet
    print(dcast(em[primary == TRUE], peak + lineages ~ statistic, value.var = "pct_matched"))
    cat("(no calls: at least two primary statistics are needed per peak)\n")
  }
}
if (length(null_scores)) {
  cal <- rbindlist(null_scores)
  fwrite(cal, file.path(dir_ev, sprintf("calibration%s.tsv", OUT_SUFFIX)), sep = "\t")
  cat("\nCalibration: how often the rule fires on non-peak regions\n")
  print(cal)
}

## Spearman correlations among statistics across the pseudo-peaks: how much independent
## information the "intersection of approaches" actually carries.
wide <- dcast(stats[region_type == "pseudo"], region_id ~ statistic + lineage, value.var = "value")
num  <- wide[, -1, with = FALSE]
keep <- names(num)[colSums(!is.na(num)) > 500]
if (length(keep) > 1) {
  cm <- cor(as.matrix(num[, ..keep]), method = "spearman", use = "pairwise.complete.obs")
  fwrite(data.table(statistic = rownames(cm), as.data.table(round(cm, 3))),
         file.path(dir_ev, sprintf("statistic_correlations%s.tsv", OUT_SUFFIX)), sep = "\t")
  cat("\nSpearman correlations among statistics across pseudo-peaks:\n")
  print(round(cm, 2))
}

cat("\nWrote evidence_matrix.tsv, calibration.tsv and statistic_correlations.tsv to\n  ",
    dir_ev, "\n", sep = "")
