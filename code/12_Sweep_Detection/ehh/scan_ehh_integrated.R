#!/usr/bin/env Rscript
## scan_ehh_integrated.R
## Compute genome-wide inES (per pop) and Rsb (between pops) from paired VCF
##
## This is the genome-wide equivalent of analyze_peak_pair_integrated.R
##
## Usage:
##   Rscript scan_ehh_integrated.R \
##     --pair BP1_TP1 \
##     --chromosome chr6 \
##     --vcf_file path/to/BP1_TP1.polarized.chr6.vcf.gz \
##     --bp_list path/to/BP1.txt \
##     --tp_list path/to/TP1.txt \
##     --output_dir path/to/output \
##     --min_maf 0.05
##
## Author: Jason Gallant Lab
## Date: 2024

suppressPackageStartupMessages({
  library(optparse)
  library(rehh)
  library(tidyverse)
})

# ========================================
# Parse Arguments
# ========================================

option_list <- list(
  make_option("--pair", type="character",
              help="Pair name (e.g., BP1_TP1, BP2_TP2, BP3_TP2)"),
  make_option("--chromosome", type="character",
              help="Chromosome name (e.g., chr6, chr8)"),
  make_option("--vcf_file", type="character",
              help="Path to paired polarized VCF file (per-chromosome)"),
  make_option("--bp_list", type="character",
              help="Path to BP population ID list"),
  make_option("--tp_list", type="character",
              help="Path to TP population ID list"),
  make_option("--output_dir", type="character",
              help="Output directory for results"),
  make_option("--min_maf", type="numeric", default=0.05,
              help="Minimum minor allele frequency filter [default: 0.05]")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Validate required arguments
required <- c("pair", "chromosome", "vcf_file", "bp_list", "tp_list", "output_dir")
missing <- required[sapply(required, function(x) is.null(opt[[x]]))]
if (length(missing) > 0) {
  stop("ERROR: Missing required arguments: ", paste(missing, collapse=", "))
}

cat("========================================\n")
cat("Integrated EHH Scan (Paired Populations)\n")
cat("========================================\n")
cat("Pair:       ", opt$pair, "\n")
cat("Chromosome: ", opt$chromosome, "\n")
cat("VCF file:   ", opt$vcf_file, "\n")
cat("BP list:    ", opt$bp_list, "\n")
cat("TP list:    ", opt$tp_list, "\n")
cat("Output dir: ", opt$output_dir, "\n")
cat("Min MAF:    ", opt$min_maf, "\n")
cat("========================================\n\n")

# ========================================
# Helper Functions (from analyze_peak_pair_integrated.R)
# ========================================

.find_col <- function(nms, candidates) {
  hit <- dplyr::first(intersect(candidates, nms))
  if (is.na(hit)) stop("Could not find any of: ", paste(candidates, collapse=", "))
  hit
}

.ind_from_hap <- function(haps) sub("(\\.|_)?[12]$", "", haps)

read_id_list <- function(path_or_vec) {
  if (length(path_or_vec) == 1 && file.exists(path_or_vec)) {
    readr::read_table(path_or_vec, col_names = FALSE, show_col_types = FALSE)[[1]]
  } else {
    as.character(path_or_vec)
  }
}

split_hh_by_groups <- function(hh_all, bp_ids, tp_ids, keep_unassigned = FALSE) {
  # normalize inputs
  bp_ids <- read_id_list(bp_ids)
  tp_ids <- read_id_list(tp_ids)

  haps <- hap.names(hh_all)
  indiv <- .ind_from_hap(haps)

  # map each haplotype to BP/TP/other
  group <- dplyr::case_when(
    indiv %in% bp_ids ~ "BP",
    indiv %in% tp_ids ~ "TP",
    TRUE              ~ "OTHER"
  )

  if (!keep_unassigned && any(group == "OTHER")) {
    message(sum(group == "OTHER"), " haplotypes are unassigned and will be ignored.")
  }

  # haplotype indices for each group
  bp_idx <- which(group == "BP")
  tp_idx <- which(group == "TP")

  # build subsets on the SAME marker set
  hh_BP <- subset(hh_all, select.hap = bp_idx)
  hh_TP <- subset(hh_all, select.hap = tp_idx)

  # tidy map (named data frame) you can reuse anywhere
  map_df <- tibble::tibble(
    hap_index = seq_along(haps),
    hap_name  = haps,
    indiv     = indiv,
    group     = group
  )

  list(
    hh   = list(BP = hh_BP, TP = hh_TP),
    idx  = list(BP = bp_idx, TP = tp_idx),
    map  = map_df
  )
}

# ========================================
# Main Analysis Pipeline
# ========================================

# 1. Load paired VCF
cat("Step 1: Loading paired VCF...\n")
if (!file.exists(opt$vcf_file)) {
  stop("ERROR: VCF file not found: ", opt$vcf_file)
}

hh_full <- tryCatch({
  data2haplohh(
    hap_file = opt$vcf_file,
    polarize_vcf = TRUE,
    chr.name = opt$chromosome,
    vcf_reader = "data.table"
  )
}, error = function(e) {
  cat("ERROR loading VCF:\n")
  cat("  ", conditionMessage(e), "\n")
  stop("Failed to load haplohh from VCF")
})

cat("  Loaded", nhap(hh_full), "haplotypes\n")
cat("  Initial markers:", nmrk(hh_full), "\n\n")

# 2. Apply MAF filter to full dataset
cat("Step 2: Filtering to MAF >=", opt$min_maf, "...\n")
hh_full <- subset(hh_full, min_maf = opt$min_maf)
cat("  Filtered markers:", nmrk(hh_full), "\n")

if (nmrk(hh_full) == 0) {
  stop("ERROR: No markers remaining after MAF filter!")
}
cat("\n")

# 3. Split into BP and TP subsets on SAME marker set
cat("Step 3: Splitting into BP and TP subsets...\n")
if (!file.exists(opt$bp_list)) {
  stop("ERROR: BP list file not found: ", opt$bp_list)
}
if (!file.exists(opt$tp_list)) {
  stop("ERROR: TP list file not found: ", opt$tp_list)
}

split_result <- split_hh_by_groups(hh_full, opt$bp_list, opt$tp_list)
hh_BP <- split_result$hh$BP
hh_TP <- split_result$hh$TP

cat("  BP haplotypes:", nhap(hh_BP), "\n")
cat("  TP haplotypes:", nhap(hh_TP), "\n")
cat("  Shared markers:", nmrk(hh_BP), "(BP) =", nmrk(hh_TP), "(TP)\n")

if (nmrk(hh_BP) != nmrk(hh_TP)) {
  stop("ERROR: Marker count mismatch after split!")
}
cat("\n")

# 4. Run scan_hh on both subsets
cat("Step 4: Running scan_hh() for BP population...\n")
scan_BP <- tryCatch({
  scan_hh(hh_BP, polarized = TRUE, discard_integration_at_border = TRUE)
}, error = function(e) {
  cat("ERROR during scan_hh (BP):\n")
  cat("  ", conditionMessage(e), "\n")
  stop("Failed to run scan_hh() for BP")
})
cat("  Computed inES for", nrow(scan_BP), "positions\n\n")

cat("Step 5: Running scan_hh() for TP population...\n")
scan_TP <- tryCatch({
  scan_hh(hh_TP, polarized = TRUE, discard_integration_at_border = TRUE)
}, error = function(e) {
  cat("ERROR during scan_hh (TP):\n")
  cat("  ", conditionMessage(e), "\n")
  stop("Failed to run scan_hh() for TP")
})
cat("  Computed inES for", nrow(scan_TP), "positions\n\n")

# 5. Extract population labels from pair name (e.g., "BP1_TP1" -> "BP1", "TP1")
pair_parts <- strsplit(opt$pair, "_")[[1]]
if (length(pair_parts) != 2) {
  stop("ERROR: Pair name must be in format BP_TP (e.g., BP1_TP1)")
}
label_BP <- pair_parts[1]
label_TP <- pair_parts[2]

cat("Step 6: Computing Rsb from scan results...\n")
rsb_result <- tryCatch({
  ines2rsb(
    scan_pop1 = scan_BP,
    scan_pop2 = scan_TP,
    popname1 = label_BP,
    popname2 = label_TP
  )
}, error = function(e) {
  cat("ERROR during ines2rsb:\n")
  cat("  ", conditionMessage(e), "\n")
  stop("Failed to compute Rsb")
})
cat("  Computed Rsb for", nrow(rsb_result), "positions\n\n")

# ========================================
# Format and Save Results
# ========================================

cat("Step 7: Formatting output data frames...\n")

# BP inES
scan_BP_df <- tibble::as_tibble(scan_BP)
chr_col <- .find_col(names(scan_BP_df), c("CHR","chr","CHROM","chrom"))
pos_col <- .find_col(names(scan_BP_df), c("POSITION","position","POS","pos"))
ines_col <- .find_col(names(scan_BP_df), c("inES","INES","ines"))
ihs_col <- dplyr::first(intersect(c("iHS","IHS","ihs"), names(scan_BP_df)))

# Standardize iHS if not computed
if (is.na(ihs_col)) {
  scan_BP_df <- scan_BP_df %>%
    dplyr::mutate(iHS = (!!sym(ines_col) - mean(!!sym(ines_col), na.rm=TRUE)) /
                        sd(!!sym(ines_col), na.rm=TRUE))
  ihs_col <- "iHS"
}

output_BP <- scan_BP_df %>%
  dplyr::select(CHR = !!sym(chr_col),
                POSITION = !!sym(pos_col),
                inES = !!sym(ines_col),
                iHS = !!sym(ihs_col))

# TP inES
scan_TP_df <- tibble::as_tibble(scan_TP)
chr_col <- .find_col(names(scan_TP_df), c("CHR","chr","CHROM","chrom"))
pos_col <- .find_col(names(scan_TP_df), c("POSITION","position","POS","pos"))
ines_col <- .find_col(names(scan_TP_df), c("inES","INES","ines"))
ihs_col <- dplyr::first(intersect(c("iHS","IHS","ihs"), names(scan_TP_df)))

if (is.na(ihs_col)) {
  scan_TP_df <- scan_TP_df %>%
    dplyr::mutate(iHS = (!!sym(ines_col) - mean(!!sym(ines_col), na.rm=TRUE)) /
                        sd(!!sym(ines_col), na.rm=TRUE))
  ihs_col <- "iHS"
}

output_TP <- scan_TP_df %>%
  dplyr::select(CHR = !!sym(chr_col),
                POSITION = !!sym(pos_col),
                inES = !!sym(ines_col),
                iHS = !!sym(ihs_col))

# Rsb
rsb_df <- tibble::as_tibble(rsb_result)
chr_col <- .find_col(names(rsb_df), c("CHR","chr","CHROM","chrom"))
pos_col <- .find_col(names(rsb_df), c("POSITION","position","POS","pos"))
rsb_col <- dplyr::first(grep("^RSB", names(rsb_df), value = TRUE))

if (is.na(rsb_col)) {
  stop("ERROR: Could not find RSB column in rsb output")
}

output_RSB <- rsb_df %>%
  dplyr::select(CHR = !!sym(chr_col),
                POSITION = !!sym(pos_col),
                RSB = !!sym(rsb_col))

# ========================================
# Save Output Files
# ========================================

cat("\nStep 8: Saving output files...\n")

# Create output directory
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

# Define output file names
# Format: {PAIR}.BP.{CHR}.ines.txt.gz  (BP population inES)
#         {PAIR}.TP.{CHR}.ines.txt.gz  (TP population inES)
#         {PAIR}.{CHR}.rsb.txt.gz       (Rsb between populations)

file_BP <- file.path(opt$output_dir, paste0(opt$pair, ".BP.", opt$chromosome, ".ines.txt.gz"))
file_TP <- file.path(opt$output_dir, paste0(opt$pair, ".TP.", opt$chromosome, ".ines.txt.gz"))
file_RSB <- file.path(opt$output_dir, paste0(opt$pair, ".", opt$chromosome, ".rsb.txt.gz"))

# Save files
readr::write_tsv(output_BP, file_BP)
readr::write_tsv(output_TP, file_TP)
readr::write_tsv(output_RSB, file_RSB)

cat("  Saved BP inES:  ", file_BP, "\n")
cat("    Size:", round(file.size(file_BP) / 1024^2, 2), "MB,",
    nrow(output_BP), "positions\n")

cat("  Saved TP inES:  ", file_TP, "\n")
cat("    Size:", round(file.size(file_TP) / 1024^2, 2), "MB,",
    nrow(output_TP), "positions\n")

cat("  Saved Rsb:      ", file_RSB, "\n")
cat("    Size:", round(file.size(file_RSB) / 1024^2, 2), "MB,",
    nrow(output_RSB), "positions\n")

# ========================================
# Summary Statistics
# ========================================

cat("\n========================================\n")
cat("Summary Statistics\n")
cat("========================================\n")

cat("\nBP Population (", label_BP, "):\n", sep="")
cat("  inES: mean =", round(mean(output_BP$inES, na.rm=TRUE), 3),
    ", SD =", round(sd(output_BP$inES, na.rm=TRUE), 3), "\n")
cat("  iHS:  mean =", round(mean(output_BP$iHS, na.rm=TRUE), 3),
    "(~0), SD =", round(sd(output_BP$iHS, na.rm=TRUE), 3), "(~1)\n")

cat("\nTP Population (", label_TP, "):\n", sep="")
cat("  inES: mean =", round(mean(output_TP$inES, na.rm=TRUE), 3),
    ", SD =", round(sd(output_TP$inES, na.rm=TRUE), 3), "\n")
cat("  iHS:  mean =", round(mean(output_TP$iHS, na.rm=TRUE), 3),
    "(~0), SD =", round(sd(output_TP$iHS, na.rm=TRUE), 3), "(~1)\n")

cat("\nRsb (", label_BP, " / ", label_TP, "):\n", sep="")
cat("  Mean:", round(mean(output_RSB$RSB, na.rm=TRUE), 3), "(should be ~0)\n")
cat("  SD:  ", round(sd(output_RSB$RSB, na.rm=TRUE), 3), "\n")
cat("  Range: [", round(min(output_RSB$RSB, na.rm=TRUE), 2), ", ",
    round(max(output_RSB$RSB, na.rm=TRUE), 2), "]\n", sep="")
cat("  |Rsb| > 2:", sum(abs(output_RSB$RSB) > 2, na.rm=TRUE), "positions\n")
cat("  |Rsb| > 3:", sum(abs(output_RSB$RSB) > 3, na.rm=TRUE), "positions\n")

cat("\n========================================\n")
cat("Done!\n")
cat("========================================\n")
