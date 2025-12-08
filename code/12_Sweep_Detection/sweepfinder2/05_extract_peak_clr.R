#!/usr/bin/env Rscript
## 05_extract_peak_clr.R
## Extract CLR values within GWAS peak regions
##
## This script reads genome-wide CLR distributions and extracts
## CLR values within the 6 major GWAS peaks for all populations.
##
## Usage:
##   Rscript code/12_Sweep_Detection/05_extract_peak_clr.R
##
## Prerequisites:
##   - 04_merge_clr_results.sh completed
##
## Author: Jason Gallant Lab
## Date: 2024

library(tidyverse)

# Find project root
root <- rprojroot::find_root(".git/index")

cat("======================================== \n")
cat("Extracting GWAS Peak CLR Values\n")
cat("======================================== \n\n")

# Define populations
pops <- c("BP1", "TP1", "BP2", "TP2", "BP3")

# Load GWAS peaks
cat("Loading GWAS peaks...\n")
peaks_file <- file.path(root, "output_data/06_Association/PKINGS_ALL_WOB_EXCLUDED/PKINGS_ALL_WOB_EXCLUDED_PEAKS.bed")

if (!file.exists(peaks_file)) {
  stop("ERROR: GWAS peaks BED file not found: ", peaks_file)
}

peaks <- read.table(peaks_file, header = FALSE,
                    col.names = c("chromosome", "start", "end", "peak", "max_log_p", "idx_bp"))

# Add chr prefix if needed
if (!grepl("^chr", peaks$chromosome[1])) {
  peaks$chromosome <- paste0("chr", peaks$chromosome)
}

cat("  Found", nrow(peaks), "GWAS peaks\n")
cat("  Chromosomes:", paste(unique(peaks$chromosome), collapse = ", "), "\n\n")

# Initialize results list
all_peak_clr <- list()

# Load genome-wide CLR for each population
for (pop in pops) {
  cat("========================================\n")
  cat("Processing population:", pop, "\n")
  cat("========================================\n")

  clr_file <- file.path(root, "output_data/12_Sweep_Detection/sweepfinder2/merged",
                        paste0(pop, ".genome_wide_clr.txt"))

  if (!file.exists(clr_file)) {
    warning("CLR file not found for ", pop, ": ", clr_file)
    next
  }

  cat("Loading:", clr_file, "\n")

  # Read genome-wide CLR
  clr_data <- read_tsv(clr_file, show_col_types = FALSE)

  cat("  Total positions:", nrow(clr_data), "\n")
  cat("  Chromosomes:", length(unique(clr_data$chromosome)), "\n")

  # Extract CLR values within each peak
  peak_clr_list <- list()

  for (i in 1:nrow(peaks)) {
    peak_chr <- peaks$chromosome[i]
    peak_start <- peaks$start[i]
    peak_end <- peaks$end[i]
    peak_num <- peaks$peak[i]

    cat("  Extracting peak", peak_num, "(", peak_chr, ":", peak_start, "-", peak_end, ")... ")

    # Filter CLR data to peak region
    peak_clr <- clr_data %>%
      filter(chromosome == peak_chr,
             position >= peak_start,
             position <= peak_end) %>%
      mutate(peak = peak_num,
             population = pop)

    cat(nrow(peak_clr), "positions\n")

    peak_clr_list[[i]] <- peak_clr
  }

  # Combine all peaks for this population
  pop_peak_clr <- bind_rows(peak_clr_list)
  all_peak_clr[[pop]] <- pop_peak_clr

  cat("  Total peak positions:", nrow(pop_peak_clr), "\n\n")
}

# Combine all populations
cat("========================================\n")
cat("Combining results\n")
cat("========================================\n")

combined_peak_clr <- bind_rows(all_peak_clr)

cat("Total rows in combined data:", nrow(combined_peak_clr), "\n")
cat("Populations:", paste(unique(combined_peak_clr$population), collapse = ", "), "\n")
cat("Peaks:", paste(unique(combined_peak_clr$peak), collapse = ", "), "\n\n")

# Save combined peak CLR values
outdir <- file.path(root, "output_data/12_Sweep_Detection/sweepfinder2")
outfile <- file.path(outdir, "peak_clr_values.txt")
write_tsv(combined_peak_clr, outfile)
cat("Saved peak CLR values:", outfile, "\n")
cat("  File size:", file.size(outfile) / 1024^2, "MB\n\n")

# Calculate summary statistics per peak per population
cat("========================================\n")
cat("Calculating summary statistics\n")
cat("========================================\n")

peak_summary <- combined_peak_clr %>%
  group_by(population, peak, chromosome) %>%
  summarize(
    n_positions = n(),
    mean_CLR = mean(CLR, na.rm = TRUE),
    median_CLR = median(CLR, na.rm = TRUE),
    max_CLR = max(CLR, na.rm = TRUE),
    sd_CLR = sd(CLR, na.rm = TRUE),
    position_max = position[which.max(CLR)],
    .groups = "drop"
  ) %>%
  arrange(desc(max_CLR))

# Save summary
outfile_summary <- file.path(outdir, "peak_clr_summary.txt")
write_tsv(peak_summary, outfile_summary)
cat("Saved summary statistics:", outfile_summary, "\n\n")

# Print summary table
cat("Peak CLR Summary (top 20 by max CLR):\n")
print(peak_summary %>% head(20), n = 20)

cat("\n")
cat("========================================\n")
cat("Peak CLR extraction complete\n")
cat("========================================\n")
cat("\nNext step: Render analysis notebook with sweepfinder2_analysis.Rmd\n")
