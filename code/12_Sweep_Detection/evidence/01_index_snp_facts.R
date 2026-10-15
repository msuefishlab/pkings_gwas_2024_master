#!/usr/bin/env Rscript
## 01_index_snp_facts.R
## Step 1 of the sweep-evidence analysis (see DESIGN.md §3).
##
## For every GWAS peak index SNP and every lineage, record:
##   - the frequency of the BP-associated allele,
##   - whether that allele is derived or ancestral (INFO/AA from the P. szaboi outgroup),
##   - whether each statistic is applicable in that lineage.
##
## The BP-associated allele is defined from the allele frequencies themselves: at each peak it is
## the allele at higher frequency in the lineage(s) in which the BP phenotype arose than in the
## matched TP reference. Both alleles' frequencies are written out, so the choice is checkable.
##
## Applicability (DESIGN.md §3): iHS and nSL need an intermediate-frequency allele and are declared
## applicable only where the BP allele frequency is <= 0.85 in that lineage. Rsb, XP-EHH, CLR, the
## pi ratio and PBS apply at any frequency.
##
## Usage (inside gwas_tools.sif, from the repo root):
##   Rscript code/12_Sweep_Detection/evidence/01_index_snp_facts.R
##
## Input:  output_data/12_Sweep_Detection/vcfs/by_group/{BP1,BP2,BP3,TP1,TP2}.polarized.vcf.gz
##         input_data/10_Introgression/peaks.bed
## Output: output_data/12_Sweep_Detection/evidence/index_snp_facts.tsv
##         output_data/12_Sweep_Detection/evidence/applicability_map.tsv

suppressPackageStartupMessages({
  library(data.table)
  library(rprojroot)
})

## .git is a directory in a clone and a file in a linked worktree, and rprojroot's
## has_file() criterion matches only the latter. Accept either, so the script runs
## wherever the repository is checked out.
root <- find_root(has_file(".git") | has_dir(".git"))
MAX_FREQ_FOR_IHS <- 0.85   # DESIGN.md section 3

vcf_dir  <- file.path(root, "output_data/12_Sweep_Detection/vcfs/by_group")
out_dir  <- file.path(root, "output_data/12_Sweep_Detection/evidence")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lineages <- c("BP1", "TP1", "BP2", "BP3", "TP2")
# lineage in which each peak's BP allele arose, and its matched TP reference (DESIGN.md section 2)
own <- data.table(
  peak = c(1L, 2L, 2L, 3L, 3L, 4L, 5L, 6L),
  bp   = c("BP1", "BP2", "BP3", "BP2", "BP3", "BP1", "BP1", "BP2"),
  tp   = c("TP1", "TP2", "TP2", "TP2", "TP2", "TP1", "TP1", "TP2")
)

peaks <- fread(file.path(root, "input_data/10_Introgression/peaks.bed"),
               col.names = c("chrom_num", "start", "end", "peak", "score", "index_snp"))
peaks[, chrom := paste0("chr", chrom_num)]
setorder(peaks, peak)
cat("Peaks:\n"); print(peaks[, .(peak, chrom, start, end, index_snp)])

## ---------------------------------------------------------------- allele counts at the index SNPs
query_one <- function(lineage) {
  vcf <- file.path(vcf_dir, sprintf("%s.polarized.vcf.gz", lineage))
  if (!file.exists(vcf)) stop("missing VCF: ", vcf)
  regions <- peaks[, paste0(chrom, ":", index_snp, "-", index_snp)]
  cmd <- sprintf("bcftools query -r '%s' -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\t%%INFO/AA\\t%%INFO/AC\\t%%INFO/AN\\n' %s",
                 paste(regions, collapse = ","), shQuote(vcf))
  x <- fread(cmd = cmd, header = FALSE, na.strings = c(".", ""),
             col.names = c("chrom", "pos", "ref", "alt", "aa", "ac", "an"))
  if (nrow(x) == 0) stop("no index SNPs found in ", vcf)
  x[, lineage := lineage]
  x[]
}

counts <- rbindlist(lapply(lineages, query_one))
counts <- merge(counts, peaks[, .(peak, chrom, index_snp)],
                by.x = c("chrom", "pos"), by.y = c("chrom", "index_snp"), all.x = TRUE)
counts[, alt_freq := ac / an]
counts[, n_hap := an]

missing <- peaks[!peak %in% counts$peak]
if (nrow(missing)) cat("WARNING: index SNP absent from the group VCFs for peak(s): ",
                       paste(missing$peak, collapse = ", "), "\n", sep = "")

## ---------------------------------------------------------------------- the BP-associated allele
## Defined per peak: the allele at higher frequency in the own lineage(s) than in the matched TP.
## Where two lineages share a peak (P2, P3) the mean difference over them is used, so one allele is
## named per peak and the per-lineage frequencies below are all of that same allele.
freq <- dcast(counts, peak + chrom + pos + ref + alt + aa ~ lineage, value.var = "alt_freq")
diffs <- own[, .(d = mean(unlist(Map(function(b, t)
  freq[peak == .BY$peak][[b]] - freq[peak == .BY$peak][[t]], bp, tp)))), by = peak]
freq <- merge(freq, diffs, by = "peak")
freq[, bp_allele := ifelse(d >= 0, alt, ref)]
freq[, bp_is_alt := bp_allele == alt]
freq[, bp_allele_is_derived := fifelse(is.na(aa), NA, bp_allele != aa)]
freq[, ancestral_allele := aa]

facts <- melt(counts[, c("peak", "lineage", "alt_freq", "n_hap"), with = FALSE],
              id.vars = c("peak", "lineage", "n_hap"), value.name = "alt_freq")[, variable := NULL][]
facts <- merge(facts, freq[, .(peak, chrom, pos, ref, alt, ancestral_allele, bp_allele,
                               bp_is_alt, bp_allele_is_derived)], by = "peak")
facts[, bp_allele_freq := fifelse(bp_is_alt, alt_freq, 1 - alt_freq)]
facts <- merge(facts, own[, .(peak, bp)], by.x = c("peak", "lineage"), by.y = c("peak", "bp"),
               all.x = TRUE)[, role := fifelse(!is.na(lineage) & paste(peak, lineage) %in%
                                                 own[, paste(peak, bp)], "own",
                                               fifelse(paste(peak, lineage) %in% own[, paste(peak, tp)],
                                                       "matched_TP", "other"))]
facts[, ihs_applicable := bp_allele_freq <= MAX_FREQ_FOR_IHS]
setcolorder(facts, c("peak", "chrom", "pos", "lineage", "role", "ref", "alt", "ancestral_allele",
                     "bp_allele", "bp_allele_is_derived", "bp_allele_freq", "n_hap", "ihs_applicable"))
setorder(facts, peak, -role, lineage)

fwrite(facts, file.path(out_dir, "index_snp_facts.tsv"), sep = "\t")

## ------------------------------------------------------------------------- the applicability map
## One row per test unit (peak x own lineage): which statistics may be used there.
appl <- merge(own, facts[role == "own", .(peak, lineage, bp_allele_freq, bp_allele_is_derived,
                                          ihs_applicable)],
              by.x = c("peak", "bp"), by.y = c("peak", "lineage"))
appl[, `:=`(
  rsb      = TRUE,          # primary, class H
  clr      = TRUE,          # primary, class S
  pi_ratio = TRUE,          # primary, class D
  xpehh    = TRUE,
  ihs      = ihs_applicable,
  nsl      = ihs_applicable,
  pbs      = TRUE
)]
appl[, ihs_applicable := NULL]
fwrite(appl, file.path(out_dir, "applicability_map.tsv"), sep = "\t")

cat("\nIndex-SNP facts (own lineages and matched TPs):\n")
print(facts[role != "other", .(peak, lineage, role, bp_allele, derived = bp_allele_is_derived,
                               freq = round(bp_allele_freq, 3), n_hap)])
cat("\nApplicability map:\n"); print(appl)
cat("\nWrote:\n  ", file.path(out_dir, "index_snp_facts.tsv"),
    "\n  ", file.path(out_dir, "applicability_map.tsv"), "\n", sep = "")
