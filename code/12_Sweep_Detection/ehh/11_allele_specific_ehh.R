#!/usr/bin/env Rscript
## 11_allele_specific_ehh.R
## Allele-specific EHH, haplotype furcation and haplotype length at each GWAS index SNP.
##
## Fig. S7 currently shows calc_ehhs(), which is EHH summed over all haplotypes in a
## population: it describes the population, not the allele. Its legend nevertheless describes
## "haplotypes carrying the BP-associated (derived) allele" against "the TP-associated
## (ancestral) allele", which is calc_ehh() and was never run. This script runs it, so the
## revised figure can show either the curve the legend promises or a legend that matches the
## curve.
##
## Note that the two differ most where the BP allele is near fixation in its lineage: there
## the population curve and the BP-allele curve nearly coincide, and the ancestral curve rests
## on very few haplotypes. index_snp_facts.tsv carries the frequencies needed to read the
## panels honestly, and the per-allele haplotype counts are written out here as well.
##
## Usage (inside rehh.sif, from the repo root):
##   Rscript code/12_Sweep_Detection/ehh/11_allele_specific_ehh.R
##
## Output: output_data/12_Sweep_Detection/ehh/allele_specific/peak{N}_{lineage}.rds
##         output_data/12_Sweep_Detection/ehh/allele_specific/summary.tsv

suppressPackageStartupMessages({
  library(data.table)
  library(rehh)
  library(rprojroot)
})

## .git is a directory in a clone and a file in a linked worktree, and rprojroot's
## has_file() criterion matches only the latter. Accept either, so the script runs
## wherever the repository is checked out.
root <- find_root(has_file(".git") | has_dir(".git"))
WIN <- 250e3   # plotting window either side of the index SNP, as in the published figure

out_dir <- file.path(root, "output_data/12_Sweep_Detection/ehh/allele_specific")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

peaks <- fread(file.path(root, "input_data/10_Introgression/peaks.bed"),
               col.names = c("chrom_num", "start", "end", "peak", "score", "index_snp"))
peaks[, chrom := paste0("chr", chrom_num)]

## peak -> the lineages whose BP allele arose there, plus the matched TP (DESIGN.md section 2)
units <- data.table(
  peak = c(1L, 2L, 2L, 3L, 3L, 4L, 5L, 6L, 6L),
  lineage = c("BP1", "BP2", "BP3", "BP2", "BP3", "BP1", "BP1", "BP2", "BP3"),
  tp = c("TP1", "TP2", "TP2", "TP2", "TP2", "TP1", "TP1", "TP2", "TP2"),
  role = c("own", "own", "own", "own", "own", "own", "own", "own", "own_completeness"))

read_ids <- function(f) sort(unique(readLines(f)))
ind_from_hap <- function(h) sub("(\\.|_)?[12]$", "", h)

summ <- list()
for (i in seq_len(nrow(units))) {
  pk <- units$peak[i]; lg <- units$lineage[i]; tp <- units$tp[i]
  pr <- peaks[peak == pk]
  vcf <- file.path(root, "output_data/12_Sweep_Detection/by_chrom",
                   sprintf("%s.polarized.%s.vcf.gz", lg, pr$chrom))
  if (!file.exists(vcf)) { cat("skip peak", pk, lg, "- missing", vcf, "\n"); next }

  cat("\nPeak ", pk, " ", lg, " (", pr$chrom, ":", pr$index_snp, ")\n", sep = "")
  hh <- data2haplohh(hap_file = vcf, polarize_vcf = TRUE, chr.name = pr$chrom,
                     vcf_reader = "data.table", verbose = FALSE)
  hh <- subset(hh, min_maf = 0.05, verbose = FALSE)

  pos <- positions(hh)
  mrk <- which.min(abs(pos - pr$index_snp))
  if (abs(pos[mrk] - pr$index_snp) > 5000) {
    cat("  nearest scored marker is ", abs(pos[mrk] - pr$index_snp),
        " bp from the index SNP; skipping\n", sep = "")
    next
  }
  cat("  marker ", mrk, " at ", pos[mrk], " (", pos[mrk] - pr$index_snp,
      " bp from the index SNP)\n", sep = "")

  ## allele-specific EHH: separate curves for the ancestral and derived haplotypes
  eh <- tryCatch(calc_ehh(hh, mrk = mrk, include_nhaplo = TRUE, discard_integration_at_border = FALSE),
                 error = function(e) { cat("  calc_ehh failed:", conditionMessage(e), "\n"); NULL })
  fu <- tryCatch(calc_furcation(hh, mrk = mrk),
                 error = function(e) { cat("  calc_furcation failed:", conditionMessage(e), "\n"); NULL })
  hl <- if (!is.null(fu)) tryCatch(calc_haplen(fu),
                                   error = function(e) { cat("  calc_haplen failed:", conditionMessage(e), "\n"); NULL })
  ## the population curve the published figure shows, for the same marker
  es <- tryCatch(calc_ehhs(hh, mrk = mrk, include_nhaplo = TRUE),
                 error = function(e) NULL)

  saveRDS(list(peak = pk, lineage = lg, tp = tp, chrom = pr$chrom,
               index_snp = pr$index_snp, marker_pos = pos[mrk], window = WIN,
               ehh = eh, ehhs = es, furcation = fu, haplen = hl),
          file.path(out_dir, sprintf("peak%d_%s.rds", pk, lg)), compress = "xz")

  n_anc <- n_der <- NA_integer_
  if (!is.null(eh) && !is.null(eh$ehh)) {
    d <- as.data.table(eh$ehh)
    nm <- names(d)
    anc <- grep("NHAPLO_A", nm, value = TRUE); der <- grep("NHAPLO_D", nm, value = TRUE)
    if (length(anc)) n_anc <- max(d[[anc[1]]], na.rm = TRUE)
    if (length(der)) n_der <- max(d[[der[1]]], na.rm = TRUE)
  }
  summ[[length(summ) + 1]] <- data.table(
    peak = pk, lineage = lg, role = units$role[i], chrom = pr$chrom,
    index_snp = pr$index_snp, marker_pos = pos[mrk],
    n_haplotypes = nhap(hh), n_ancestral = n_anc, n_derived = n_der,
    has_ehh = !is.null(eh), has_furcation = !is.null(fu))
  cat("  haplotypes: ", nhap(hh), " (ancestral ", n_anc, ", derived ", n_der, ")\n", sep = "")
}

if (length(summ)) {
  s <- rbindlist(summ)
  fwrite(s, file.path(out_dir, "summary.tsv"), sep = "\t")
  cat("\n", strrep("=", 60), "\n", sep = ""); print(s)
  cat("\nWrote ", out_dir, "\n", sep = "")
}
