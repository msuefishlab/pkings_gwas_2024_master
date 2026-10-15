#!/usr/bin/env Rscript
## 10_genome_stats.R
## Real iHS, XP-EHH and Rsb from the full scan_hh output saved by scan_ehh_integrated.R.
##
## The published analysis reported a statistic called "iHS" that is inES standardized to mean 0
## and SD 1 within each population. That is not iHS: iHS contrasts haplotypes carrying the
## derived allele with those carrying the ancestral allele at the same site, and it is
## standardized within allele-frequency bins because haplotype length depends strongly on
## frequency. ihh2ihs() below computes the real statistic. XP-EHH was never computed at all.
##
## All three are standardized genome-wide, across all 25 chromosomes at once, so a peak is
## measured against the whole genome rather than against its own chromosome.
##
## Usage (inside rehh.sif, from the repo root):
##   Rscript code/12_Sweep_Detection/ehh/10_genome_stats.R [--subsample 16]
##
## Input:  ${sweep_scratch}/12_Sweep_Detection/ehh_rds/{PAIR}.{BP,TP}[.nN].{chrom}.scan.rds
## Output: output_data/12_Sweep_Detection/ehh/{ihs,xpehh,rsb_check}/{lineage}.genome_wide.txt.gz
##         output_data/12_Sweep_Detection/ehh/genome_stats_summary.tsv

suppressPackageStartupMessages({
  library(data.table)
  library(rehh)
  library(rprojroot)
})

root <- find_root(".git")
args <- commandArgs(trailingOnly = TRUE)
subsample <- if (length(args) >= 2 && args[1] == "--subsample") as.integer(args[2]) else NA_integer_
tag <- if (is.na(subsample)) "" else paste0(".n", subsample)

rds_dir <- Sys.getenv("SWEEP_SCRATCH", "/mnt/gs21/scratch/jgallant/pkings_sweeps")
rds_dir <- file.path(rds_dir, "12_Sweep_Detection/ehh_rds")
out_base <- file.path(root, "output_data/12_Sweep_Detection/ehh")

pairs <- c(BP1 = "BP1_TP1", BP2 = "BP2_TP2", BP3 = "BP3_TP2")
tp_of <- c(BP1 = "TP1",     BP2 = "TP2",     BP3 = "TP2")
chroms <- paste0("chr", 1:25)

## Read the 25 per-chromosome scans for one group and stack them. rehh's standardizing
## functions expect a single data frame covering everything they should standardize over.
read_scan <- function(pair, grp) {
  files <- file.path(rds_dir, sprintf("%s.%s%s.%s.scan.rds", pair, grp, tag, chroms))
  have  <- file.exists(files)
  if (!any(have)) stop("no scan RDS files for ", pair, " ", grp, " in ", rds_dir)
  if (!all(have)) cat("  NOTE:", sum(!have), "of 25 chromosomes missing for", pair, grp, "\n")
  rbindlist(lapply(files[have], readRDS), use.names = TRUE, fill = TRUE)
}

write_stat <- function(dt, stat, lineage, value_col) {
  d <- as.data.table(dt)
  setnames(d, c("CHR", "POSITION"), c("chrom", "pos"), skip_absent = TRUE)
  d <- d[, .(chrom, pos, value = get(value_col))][!is.na(value)]
  dir <- file.path(out_base, stat); dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(dir, sprintf("%s%s.genome_wide.txt.gz", lineage, tag))
  fwrite(d, f, sep = "\t")
  cat(sprintf("  wrote %-58s %8d sites\n", f, nrow(d)))
  data.table(statistic = stat, lineage = lineage, subsample = subsample, n_sites = nrow(d),
             mean = mean(d$value), sd = sd(d$value),
             p95 = quantile(d$value, 0.95, na.rm = TRUE),
             p99 = quantile(d$value, 0.99, na.rm = TRUE))
}

summ <- list()
for (bp in names(pairs)) {
  pair <- pairs[[bp]]; tp <- tp_of[[bp]]
  cat("\n", strrep("=", 64), "\n", pair, if (nzchar(tag)) paste0(" (", tag, ")") else "", "\n", sep = "")

  scan_bp <- read_scan(pair, "BP")
  scan_tp <- read_scan(pair, "TP")
  cat("  BP sites:", nrow(scan_bp), " TP sites:", nrow(scan_tp), "\n")

  ## --- real iHS, standardized in allele-frequency bins, genome-wide
  for (grp in c("BP", "TP")) {
    sc  <- if (grp == "BP") scan_bp else scan_tp
    lab <- if (grp == "BP") bp else tp
    if (grp == "TP" && bp == "BP3") next          # TP2 already done with BP2
    ih <- tryCatch(ihh2ihs(sc, freqbin = 0.025, verbose = FALSE),
                   error = function(e) { cat("  iHS failed for", lab, ":", conditionMessage(e), "\n"); NULL })
    if (!is.null(ih)) {
      d <- as.data.table(ih$ihs)
      ## |iHS| is what a sweep makes large, in either allele direction
      d[, ABS_IHS := abs(IHS)]
      summ[[length(summ) + 1]] <- write_stat(d, "ihs", lab, "ABS_IHS")
    }
  }

  ## --- XP-EHH, BP versus its matched TP
  xp <- tryCatch(ies2xpehh(scan_bp, scan_tp, popname1 = bp, popname2 = tp, verbose = FALSE),
                 error = function(e) { cat("  XP-EHH failed:", conditionMessage(e), "\n"); NULL })
  if (!is.null(xp)) {
    d <- as.data.table(xp)
    col <- grep("^XPEHH", names(d), value = TRUE)[1]
    summ[[length(summ) + 1]] <- write_stat(d, "xpehh", bp, col)
  }

  ## --- Rsb recomputed here, as a check that the RDS files reproduce the published scan
  rs <- tryCatch(ines2rsb(scan_bp, scan_tp, popname1 = bp, popname2 = tp, verbose = FALSE),
                 error = function(e) { cat("  Rsb failed:", conditionMessage(e), "\n"); NULL })
  if (!is.null(rs)) {
    d <- as.data.table(rs)
    col <- grep("^RSB", names(d), value = TRUE)[1]
    summ[[length(summ) + 1]] <- write_stat(d, "rsb_check", bp, col)
  }
}

if (length(summ)) {
  s <- rbindlist(summ)
  f <- file.path(out_base, sprintf("genome_stats_summary%s.tsv", tag))
  fwrite(s, f, sep = "\t")
  cat("\n", strrep("=", 64), "\nSummary:\n", sep = ""); print(s)
  cat("\nWrote ", f, "\n", sep = "")
}
