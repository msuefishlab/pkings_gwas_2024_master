#!/usr/bin/env Rscript
## scan_rsb_genome_wide.R
## Compute genome-wide Rsb from pre-computed inES values
##
## This script loads two inES files (from scan_ines_genome_wide.R),
## computes Rsb using rehh::ines2rsb(), and saves the results.
##
## Usage:
##   Rscript scan_rsb_genome_wide.R --pop1 BP1 --pop2 TP1 --chromosome chr1 \
##     --ines1 BP1.chr1.ines.txt.gz --ines2 TP1.chr1.ines.txt.gz \
##     --output_file out.txt.gz
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
  make_option("--pop1", type="character",
              help="First population name (e.g., BP1)"),
  make_option("--pop2", type="character",
              help="Second population name (e.g., TP1)"),
  make_option("--chromosome", type="character",
              help="Chromosome name (e.g., chr1, chr2)"),
  make_option("--ines1", type="character",
              help="Path to inES file for population 1"),
  make_option("--ines2", type="character",
              help="Path to inES file for population 2"),
  make_option("--output_file", type="character",
              help="Path to output file (will be gzipped)")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Validate arguments
if (is.null(opt$pop1) || is.null(opt$pop2) || is.null(opt$chromosome) ||
    is.null(opt$ines1) || is.null(opt$ines2) || is.null(opt$output_file)) {
  stop("ERROR: All arguments are required")
}

cat("========================================\n")
cat("Rsb Genome-Wide Scan (R script)\n")
cat("========================================\n")
cat("Population 1:", opt$pop1, "\n")
cat("Population 2:", opt$pop2, "\n")
cat("Chromosome:", opt$chromosome, "\n")
cat("inES file 1:", opt$ines1, "\n")
cat("inES file 2:", opt$ines2, "\n")
cat("Output file:", opt$output_file, "\n")
cat("========================================\n\n")

# Check files exist
if (!file.exists(opt$ines1)) {
  stop("ERROR: inES file not found for ", opt$pop1, ": ", opt$ines1)
}
if (!file.exists(opt$ines2)) {
  stop("ERROR: inES file not found for ", opt$pop2, ": ", opt$ines2)
}

# Load inES files
cat("Loading inES data...\n")
ines1 <- read_tsv(opt$ines1, show_col_types = FALSE)
ines2 <- read_tsv(opt$ines2, show_col_types = FALSE)

cat("  Pop1 (", opt$pop1, "):", nrow(ines1), "positions\n")
cat("  Pop2 (", opt$pop2, "):", nrow(ines2), "positions\n")

# Check that both files have same positions
# If not, we'll need to merge and match positions
if (!all(ines1$POSITION == ines2$POSITION)) {
  cat("  WARNING: Position lists differ between populations\n")
  cat("  Performing inner join on positions...\n")

  # Merge on chromosome and position
  merged <- inner_join(
    ines1 %>% select(CHR, POSITION, inES),
    ines2 %>% select(CHR, POSITION, inES),
    by = c("CHR", "POSITION"),
    suffix = c("_1", "_2")
  )

  cat("  Matched positions:", nrow(merged), "\n")

  # Reconstruct inES data frames with matched positions
  ines1 <- merged %>% select(CHR, POSITION, inES = inES_1)
  ines2 <- merged %>% select(CHR, POSITION, inES = inES_2)
}

# Convert to scan_hh format for ines2rsb()
# The function expects objects with specific structure
# We'll create minimal scan_hh-like objects

scan1 <- list(CHR = ines1$CHR,
              POSITION = ines1$POSITION,
              inES = ines1$inES)
class(scan1) <- "scanhs"  # Mark as scan object

scan2 <- list(CHR = ines2$CHR,
              POSITION = ines2$POSITION,
              inES = ines2$inES)
class(scan2) <- "scanhs"

# Compute Rsb
cat("\nComputing Rsb...\n")
rsb_result <- tryCatch({
  ines2rsb(scan_pop1 = scan1,
           scan_pop2 = scan2,
           popname1 = opt$pop1,
           popname2 = opt$pop2)
}, error = function(e) {
  cat("ERROR during ines2rsb:\n")
  cat("  ", conditionMessage(e), "\n")
  stop("Failed to compute Rsb")
})

# Convert to tibble
rsb_df <- as_tibble(rsb_result)

# Find column names (robust to different versions)
chr_col <- first(intersect(c("CHR","chr","CHROM","chrom"), names(rsb_df)))
pos_col <- first(intersect(c("POSITION","position","POS","pos"), names(rsb_df)))

# Find Rsb column (should be RSB_{POP1}_{POP2})
rsb_col <- first(grep("^RSB", names(rsb_df), value = TRUE))

if (is.na(chr_col) || is.na(pos_col) || is.na(rsb_col)) {
  cat("ERROR: Could not find required columns in Rsb output\n")
  cat("  Available columns:", paste(names(rsb_df), collapse=", "), "\n")
  stop("Column name mismatch")
}

# Create output data frame
output_df <- rsb_df %>%
  select(CHR = !!sym(chr_col),
         POSITION = !!sym(pos_col),
         !!rsb_col)  # Keep original Rsb column name

# Summary statistics
cat("\nSummary statistics:\n")
cat("  Positions:", nrow(output_df), "\n")
cat("  Rsb column:", rsb_col, "\n")
cat("  Rsb range: [", min(output_df[[rsb_col]], na.rm=TRUE), ", ",
    max(output_df[[rsb_col]], na.rm=TRUE), "]\n", sep="")
cat("  Rsb mean:", mean(output_df[[rsb_col]], na.rm=TRUE),
    "(should be ~0)\n")
cat("  Rsb SD:", sd(output_df[[rsb_col]], na.rm=TRUE), "\n")
cat("  |Rsb| > 2:", sum(abs(output_df[[rsb_col]]) > 2, na.rm=TRUE),
    "positions\n")
cat("  |Rsb| > 3:", sum(abs(output_df[[rsb_col]]) > 3, na.rm=TRUE),
    "positions\n")

# Save to compressed TSV
cat("\nSaving output to:", opt$output_file, "\n")
write_tsv(output_df, opt$output_file)

# Report file size
file_size_mb <- file.size(opt$output_file) / 1024^2
cat("  File size:", round(file_size_mb, 2), "MB\n")

cat("\nDone!\n")
cat("========================================\n")
