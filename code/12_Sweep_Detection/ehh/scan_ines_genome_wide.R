#!/usr/bin/env Rscript
## scan_ines_genome_wide.R
## Compute genome-wide inES/iHS for one population, one chromosome
##
## This script loads a polarized VCF, filters to a single chromosome,
## runs rehh::scan_hh() to compute inES values, standardizes to iHS,
## and saves the results as a compressed TSV.
##
## Usage:
##   Rscript scan_ines_genome_wide.R --population BP1 --chromosome chr1 \
##     --vcf_file path/to/vcf.gz --output_file out.txt.gz [--min_maf 0.05]
##
## Author: Jason Gallant Lab
## Date: 2024

suppressPackageStartupMessages({
  library(optparse)
  library(rehh)
  library(tidyverse)
})

# Parse command line arguments
option_list <- list(
  make_option("--population", type="character",
              help="Population name (e.g., BP1, TP1)"),
  make_option("--chromosome", type="character",
              help="Chromosome name (e.g., chr1, chr2)"),
  make_option("--vcf_file", type="character",
              help="Path to polarized VCF file"),
  make_option("--output_file", type="character",
              help="Path to output file (will be gzipped)"),
  make_option("--min_maf", type="numeric", default=0.05,
              help="Minimum minor allele frequency filter [default: 0.05]")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Validate arguments
if (is.null(opt$population) || is.null(opt$chromosome) ||
    is.null(opt$vcf_file) || is.null(opt$output_file)) {
  stop("ERROR: All arguments (--population, --chromosome, --vcf_file, --output_file) are required")
}

cat("========================================\n")
cat("inES/iHS Genome-Wide Scan (R script)\n")
cat("========================================\n")
cat("Population:", opt$population, "\n")
cat("Chromosome:", opt$chromosome, "\n")
cat("VCF file:", opt$vcf_file, "\n")
cat("Output file:", opt$output_file, "\n")
cat("Min MAF:", opt$min_maf, "\n")
cat("========================================\n\n")

# Check VCF exists
if (!file.exists(opt$vcf_file)) {
  stop("ERROR: VCF file not found: ", opt$vcf_file)
}

# Load haplohh from VCF
cat("Loading haplohh from VCF...\n")
cat("  This may take several minutes for large chromosomes\n")

hh <- tryCatch({
  data2haplohh(
    hap_file = opt$vcf_file,
    polarize_vcf = TRUE,
    chr.name = opt$chromosome,
    vcf_reader = "data.table"  # Faster VCF parsing
  )
}, error = function(e) {
  cat("ERROR loading VCF:\n")
  cat("  ", conditionMessage(e), "\n")
  stop("Failed to load haplohh from VCF")
})

cat("  Loaded", nhap(hh), "haplotypes\n")
cat("  Initial markers:", nmrk(hh), "\n")

# Apply MAF filter
cat("\nFiltering to MAF >=", opt$min_maf, "...\n")
hh <- subset(hh, min_maf = opt$min_maf)
cat("  Filtered markers:", nmrk(hh), "\n")

if (nmrk(hh) == 0) {
  stop("ERROR: No markers remaining after MAF filter!")
}

# Run scan_hh to compute inES
cat("\nRunning scan_hh()...\n")
cat("  This may take 10-30 minutes depending on chromosome size\n")

scan_result <- tryCatch({
  scan_hh(hh, polarized = TRUE, discard_integration_at_border = TRUE)
}, error = function(e) {
  cat("ERROR during scan_hh:\n")
  cat("  ", conditionMessage(e), "\n")
  stop("Failed to run scan_hh()")
})

cat("  Computed inES for", nrow(scan_result), "positions\n")

# Extract results to data frame
scan_df <- as_tibble(scan_result)

# Find column names (robust to different rehh versions)
chr_col <- first(intersect(c("CHR","chr","CHROM","chrom"), names(scan_df)))
pos_col <- first(intersect(c("POSITION","position","POS","pos"), names(scan_df)))
ines_col <- first(intersect(c("inES","INES","ines"), names(scan_df)))
ihs_col <- first(intersect(c("iHS","IHS","ihs"), names(scan_df)))

if (is.na(chr_col) || is.na(pos_col) || is.na(ines_col)) {
  cat("ERROR: Could not find required columns in scan_hh output\n")
  cat("  Available columns:", paste(names(scan_df), collapse=", "), "\n")
  stop("Column name mismatch")
}

# If iHS not computed by scan_hh, standardize inES manually
if (is.na(ihs_col)) {
  cat("\nStandardizing inES to iHS (mean=0, SD=1)...\n")
  scan_df <- scan_df %>%
    mutate(iHS = (!!sym(ines_col) - mean(!!sym(ines_col), na.rm=TRUE)) /
                 sd(!!sym(ines_col), na.rm=TRUE))
  ihs_col <- "iHS"
}

# Create output data frame
output_df <- scan_df %>%
  select(CHR = !!sym(chr_col),
         POSITION = !!sym(pos_col),
         inES = !!sym(ines_col),
         iHS = !!sym(ihs_col))

# Summary statistics
cat("\nSummary statistics:\n")
cat("  Positions:", nrow(output_df), "\n")
cat("  inES range: [", min(output_df$inES, na.rm=TRUE), ", ",
    max(output_df$inES, na.rm=TRUE), "]\n", sep="")
cat("  inES mean:", mean(output_df$inES, na.rm=TRUE), "\n")
cat("  inES SD:", sd(output_df$inES, na.rm=TRUE), "\n")
cat("  iHS range: [", min(output_df$iHS, na.rm=TRUE), ", ",
    max(output_df$iHS, na.rm=TRUE), "]\n", sep="")
cat("  iHS mean:", mean(output_df$iHS, na.rm=TRUE),
    "(should be ~0)\n")
cat("  iHS SD:", sd(output_df$iHS, na.rm=TRUE),
    "(should be ~1)\n")

# Save to compressed TSV
cat("\nSaving output to:", opt$output_file, "\n")
write_tsv(output_df, opt$output_file)

# Report file size
file_size_mb <- file.size(opt$output_file) / 1024^2
cat("  File size:", round(file_size_mb, 2), "MB\n")

cat("\nDone!\n")
cat("========================================\n")
