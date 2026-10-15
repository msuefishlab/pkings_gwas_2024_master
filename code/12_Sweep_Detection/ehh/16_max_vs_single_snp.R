#!/usr/bin/env Rscript
## Diagnostic for DESIGN.md section 1: what the published maximum-versus-single-SNP
## comparison actually measures. Not part of the pipeline; it documents the flaw it replaces.
suppressPackageStartupMessages({ library(data.table); library(rprojroot) })
root <- find_root(".git")
source(file.path(root, "code/12_Sweep_Detection/evidence/region_stats_lib.R"))

peaks <- load_peaks(root)
peaks[, width := end - start]
fai <- fread(file.path(root, "input_data/00_Reference_Genome",
  "jordan-mic4273-mb-hirise-20xx3__06-13-2024__final_assembly.fasta.renamed.fa.fai"),
  select = 1:2, col.names = c("chrom", "length"))[chrom %in% paste0("chr", 1:25)]
BUF <- 250e3

mdir <- rsb_dir(root)
load_one <- function(f, col, absval) {
  d <- fread(cmd = sprintf("zcat < %s", shQuote(file.path(mdir, f))))
  setnames(d, c("CHR", "POSITION"), c("chrom", "pos"), skip_absent = TRUE)
  d <- d[, .(chrom, pos, value = get(col))][!is.na(value)]
  if (absval) d[, value := abs(value)]
  d
}

report <- function(label, d) {
  srt <- sort(d$value)
  pctl <- function(v) 100 * findInterval(v, srt) / length(srt)
  q99 <- quantile(d$value, 0.99, names = FALSE)
  cat("\n", strrep("=", 100), "\n", label, "  |  ", format(nrow(d), big.mark = ","),
      " scored SNPs genome-wide  |  single-SNP 99th pct = ", signif(q99, 5), "\n",
      strrep("=", 100), "\n", sep = "")
  cat(sprintf("%-5s %7s %6s %11s %9s %8s | %9s %9s %9s\n", "peak", "width", "nSNP",
              "peak max", "pub pct", "1-.99^n", "null>max", "null>q99", "med pct"))
  for (i in seq_len(nrow(peaks))) {
    p <- peaks[i]; W <- p$width
    w <- d[chrom == p$chrom & pos >= p$start & pos <= p$end]
    if (!nrow(w)) next
    n <- nrow(w); mx <- max(w$value)

    ## non-peak windows of the SAME width, tiled across every chromosome
    bins <- rbindlist(lapply(seq_len(nrow(fai)), function(j) {
      ch <- fai$chrom[j]
      s <- d[chrom == ch]
      if (!nrow(s)) return(NULL)
      b <- s[, .(wmax = max(value), nsnp = .N, wstart = min(pos)),
             by = .(bin = floor((pos - 1) / W))]
      b[, `:=`(chrom = ch, lo = bin * W + 1, hi = (bin + 1) * W)]
      b
    }))
    for (k in seq_len(nrow(peaks))) bins <- bins[!(chrom == peaks$chrom[k] &
      hi > peaks$start[k] - BUF & lo < peaks$end[k] + BUF)]
    bins <- bins[nsnp > 0]
    cat(sprintf("%-5d %7s %6d %11.4g %9.3f %8.3f | %8.1f%% %8.1f%% %9.3f\n",
                p$peak, format(W, big.mark = ","), n, mx, pctl(mx), 1 - 0.99^n,
                100 * mean(bins$wmax > mx), 100 * mean(bins$wmax > q99),
                median(pctl(bins$wmax))))
  }
}

for (lg in c("BP1", "BP2", "BP3", "TP1", "TP2"))
  report(paste0("inES  ", lg), load_one(sprintf("%s.genome_wide_ines.txt.gz", lg), "inES", FALSE))
for (pr in c("BP1_TP1", "BP2_TP2", "BP3_TP2"))
  report(paste0("|Rsb|  ", pr), load_one(sprintf("%s.genome_wide_rsb.txt.gz", pr), "RSB", TRUE))

cat("\n\nSIGN CHECK: was each peak's max |Rsb| attained at a positive value (longer haplotypes\n",
    "in BP, the direction a BP sweep predicts) or a negative one (longer in TP)?\n", sep = "")
for (pr in c("BP1_TP1", "BP2_TP2", "BP3_TP2")) {
  d <- load_one(sprintf("%s.genome_wide_rsb.txt.gz", pr), "RSB", FALSE)
  cat("\n", pr, "\n", sep = "")
  for (i in seq_len(nrow(peaks))) {
    p <- peaks[i]
    w <- d[chrom == p$chrom & pos >= p$start & pos <= p$end]
    if (!nrow(w)) next
    j <- which.max(abs(w$value))
    cat(sprintf("  peak %d: max|Rsb| = %7.4f at %s:%d  -> signed Rsb %+8.4f  (%s)\n",
                p$peak, abs(w$value[j]), p$chrom, w$pos[j], w$value[j],
                if (w$value[j] > 0) "BP longer" else "TP LONGER - wrong direction"))
  }
}
