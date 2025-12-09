#!/usr/bin/env Rscript
## 04_extract_peak_ehh.R
## Extract EHH statistics within GWAS peaks and calculate percentile rankings
##
## This script performs three main tasks:
## 1. Calculate genome-wide summary statistics (mean, SD, percentiles) for inES/iHS/Rsb
## 2. Extract EHH values within the 6 major GWAS peaks for all populations
## 3. Calculate percentile rankings and empirical p-values for peak maxima
##
## Usage:
##   Rscript code/12_Sweep_Detection/ehh/04_extract_peak_ehh.R
##
## Prerequisites:
##   - 03_merge_ehh_results.sh completed (genome-wide merged files must exist)
##
## Author: Jason Gallant Lab
## Date: 2024

suppressPackageStartupMessages({
  library(tidyverse)
  library(rprojroot)
})

# Find project root
root <- find_root(".git/index")

cat("======================================== \n")
cat("Extracting GWAS Peak EHH Values\n")
cat("======================================== \n\n")

# Define populations and comparisons
pops <- c("BP1", "TP1", "BP2", "TP2", "BP3")
comparisons <- c("BP1_TP1", "BP2_TP2", "BP3_TP2")

# Define directories
merged_dir <- file.path(root, "output_data/12_Sweep_Detection/ehh/merged")
outdir <- file.path(root, "output_data/12_Sweep_Detection/ehh")

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

# ======================================================================
# Part 1: Calculate Genome-Wide Summary Statistics for inES/iHS
# ======================================================================

cat("========================================\n")
cat("Part 1: Genome-Wide inES/iHS Statistics\n")
cat("========================================\n\n")

ines_summary_list <- list()

for (pop in pops) {
  cat("Processing population:", pop, "\n")

  ines_file <- file.path(merged_dir, paste0(pop, ".genome_wide_ines.txt.gz"))

  if (!file.exists(ines_file)) {
    warning("  inES file not found: ", ines_file, "\n")
    next
  }

  cat("  Loading:", basename(ines_file), "\n")

  # Read genome-wide inES/iHS
  ines_data <- read_tsv(ines_file, show_col_types = FALSE)

  cat("    Total positions:", nrow(ines_data), "\n")

  # Calculate summary statistics for inES
  ines_stats <- tibble(
    population = pop,
    metric = "inES",
    n_snps = nrow(ines_data),
    mean = mean(ines_data$inES, na.rm = TRUE),
    sd = sd(ines_data$inES, na.rm = TRUE),
    min = min(ines_data$inES, na.rm = TRUE),
    max = max(ines_data$inES, na.rm = TRUE),
    p01 = quantile(ines_data$inES, 0.01, na.rm = TRUE),
    p05 = quantile(ines_data$inES, 0.05, na.rm = TRUE),
    p10 = quantile(ines_data$inES, 0.10, na.rm = TRUE),
    p50 = quantile(ines_data$inES, 0.50, na.rm = TRUE),
    p90 = quantile(ines_data$inES, 0.90, na.rm = TRUE),
    p95 = quantile(ines_data$inES, 0.95, na.rm = TRUE),
    p99 = quantile(ines_data$inES, 0.99, na.rm = TRUE)
  )

  # Calculate summary statistics for iHS
  ihs_stats <- tibble(
    population = pop,
    metric = "iHS",
    n_snps = nrow(ines_data),
    mean = mean(ines_data$iHS, na.rm = TRUE),
    sd = sd(ines_data$iHS, na.rm = TRUE),
    min = min(ines_data$iHS, na.rm = TRUE),
    max = max(ines_data$iHS, na.rm = TRUE),
    p01 = quantile(ines_data$iHS, 0.01, na.rm = TRUE),
    p05 = quantile(ines_data$iHS, 0.05, na.rm = TRUE),
    p10 = quantile(ines_data$iHS, 0.10, na.rm = TRUE),
    p50 = quantile(ines_data$iHS, 0.50, na.rm = TRUE),
    p90 = quantile(ines_data$iHS, 0.90, na.rm = TRUE),
    p95 = quantile(ines_data$iHS, 0.95, na.rm = TRUE),
    p99 = quantile(ines_data$iHS, 0.99, na.rm = TRUE)
  )

  ines_summary_list[[pop]] <- bind_rows(ines_stats, ihs_stats)

  cat("    inES mean:", round(ines_stats$mean, 4),
      " SD:", round(ines_stats$sd, 4), "\n")
  cat("    iHS mean:", round(ihs_stats$mean, 4),
      " (should be ~0) SD:", round(ihs_stats$sd, 4), " (should be ~1)\n\n")
}

# Combine and save inES/iHS summary
ines_summary <- bind_rows(ines_summary_list)
ines_summary_file <- file.path(outdir, "genome_wide_ines_summary.txt")
write_tsv(ines_summary, ines_summary_file)
cat("Saved inES/iHS summary:", ines_summary_file, "\n\n")

# ======================================================================
# Part 2: Calculate Genome-Wide Summary Statistics for Rsb
# ======================================================================

cat("========================================\n")
cat("Part 2: Genome-Wide Rsb Statistics\n")
cat("========================================\n\n")

rsb_summary_list <- list()

for (comp in comparisons) {
  cat("Processing comparison:", comp, "\n")

  rsb_file <- file.path(merged_dir, paste0(comp, ".genome_wide_rsb.txt.gz"))

  if (!file.exists(rsb_file)) {
    warning("  Rsb file not found: ", rsb_file, "\n")
    next
  }

  cat("  Loading:", basename(rsb_file), "\n")

  # Read genome-wide Rsb
  rsb_data <- read_tsv(rsb_file, show_col_types = FALSE)

  # Find Rsb column (should be RSB_{POP1}_{POP2})
  rsb_col <- grep("^RSB", names(rsb_data), value = TRUE)[1]

  if (is.na(rsb_col)) {
    warning("  Could not find Rsb column in ", comp, "\n")
    next
  }

  cat("    Total positions:", nrow(rsb_data), "\n")
  cat("    Rsb column:", rsb_col, "\n")

  # Calculate summary statistics
  rsb_stats <- tibble(
    comparison = comp,
    metric = "Rsb",
    n_snps = nrow(rsb_data),
    mean = mean(rsb_data[[rsb_col]], na.rm = TRUE),
    sd = sd(rsb_data[[rsb_col]], na.rm = TRUE),
    min = min(rsb_data[[rsb_col]], na.rm = TRUE),
    max = max(rsb_data[[rsb_col]], na.rm = TRUE),
    p01 = quantile(rsb_data[[rsb_col]], 0.01, na.rm = TRUE),
    p05 = quantile(rsb_data[[rsb_col]], 0.05, na.rm = TRUE),
    p10 = quantile(rsb_data[[rsb_col]], 0.10, na.rm = TRUE),
    p50 = quantile(rsb_data[[rsb_col]], 0.50, na.rm = TRUE),
    p90 = quantile(rsb_data[[rsb_col]], 0.90, na.rm = TRUE),
    p95 = quantile(rsb_data[[rsb_col]], 0.95, na.rm = TRUE),
    p99 = quantile(rsb_data[[rsb_col]], 0.99, na.rm = TRUE)
  )

  rsb_summary_list[[comp]] <- rsb_stats

  cat("    Rsb mean:", round(rsb_stats$mean, 4),
      " (should be ~0) SD:", round(rsb_stats$sd, 4), "\n\n")
}

# Combine and save Rsb summary
rsb_summary <- bind_rows(rsb_summary_list)
rsb_summary_file <- file.path(outdir, "genome_wide_rsb_summary.txt")
write_tsv(rsb_summary, rsb_summary_file)
cat("Saved Rsb summary:", rsb_summary_file, "\n\n")

# ======================================================================
# Part 3: Extract Peak Values for inES/iHS
# ======================================================================

cat("========================================\n")
cat("Part 3: Extracting Peak inES/iHS Values\n")
cat("========================================\n\n")

all_peak_ines <- list()

for (pop in pops) {
  cat("Processing population:", pop, "\n")

  ines_file <- file.path(merged_dir, paste0(pop, ".genome_wide_ines.txt.gz"))

  if (!file.exists(ines_file)) {
    warning("  inES file not found: ", ines_file, "\n")
    next
  }

  # Read genome-wide inES/iHS
  ines_data <- read_tsv(ines_file, show_col_types = FALSE)

  # Extract values within each peak
  peak_ines_list <- list()

  for (i in 1:nrow(peaks)) {
    peak_chr <- peaks$chromosome[i]
    peak_start <- peaks$start[i]
    peak_end <- peaks$end[i]
    peak_num <- peaks$peak[i]

    cat("  Extracting peak", peak_num, "(", peak_chr, ":", peak_start, "-", peak_end, ")... ")

    # Filter inES data to peak region
    peak_ines <- ines_data %>%
      filter(CHR == peak_chr,
             POSITION >= peak_start,
             POSITION <= peak_end) %>%
      mutate(peak = peak_num,
             population = pop)

    cat(nrow(peak_ines), "positions\n")

    peak_ines_list[[i]] <- peak_ines
  }

  # Combine all peaks for this population
  pop_peak_ines <- bind_rows(peak_ines_list)
  all_peak_ines[[pop]] <- pop_peak_ines

  cat("  Total peak positions:", nrow(pop_peak_ines), "\n\n")
}

# Combine all populations
combined_peak_ines <- bind_rows(all_peak_ines)

cat("Total rows in combined inES/iHS data:", nrow(combined_peak_ines), "\n")
cat("Populations:", paste(unique(combined_peak_ines$population), collapse = ", "), "\n")
cat("Peaks:", paste(unique(combined_peak_ines$peak), collapse = ", "), "\n\n")

# Save combined peak inES values
ines_outfile <- file.path(outdir, "peak_ines_values.txt")
write_tsv(combined_peak_ines, ines_outfile)
cat("Saved peak inES values:", ines_outfile, "\n")
cat("  File size:", round(file.size(ines_outfile) / 1024^2, 2), "MB\n\n")

# Calculate summary statistics per peak per population
peak_ines_summary <- combined_peak_ines %>%
  group_by(population, peak, CHR) %>%
  summarize(
    n_positions = n(),
    mean_inES = mean(inES, na.rm = TRUE),
    median_inES = median(inES, na.rm = TRUE),
    max_inES = max(inES, na.rm = TRUE),
    sd_inES = sd(inES, na.rm = TRUE),
    position_max_inES = POSITION[which.max(inES)],
    mean_iHS = mean(iHS, na.rm = TRUE),
    median_iHS = median(iHS, na.rm = TRUE),
    max_iHS = max(iHS, na.rm = TRUE),
    min_iHS = min(iHS, na.rm = TRUE),
    max_abs_iHS = max(abs(iHS), na.rm = TRUE),
    position_max_abs_iHS = POSITION[which.max(abs(iHS))],
    .groups = "drop"
  ) %>%
  arrange(desc(max_abs_iHS))

# Save summary
ines_summary_outfile <- file.path(outdir, "peak_ines_summary.txt")
write_tsv(peak_ines_summary, ines_summary_outfile)
cat("Saved inES/iHS summary:", ines_summary_outfile, "\n\n")

# ======================================================================
# Part 4: Extract Peak Values for Rsb
# ======================================================================

cat("========================================\n")
cat("Part 4: Extracting Peak Rsb Values\n")
cat("========================================\n\n")

all_peak_rsb <- list()

for (comp in comparisons) {
  cat("Processing comparison:", comp, "\n")

  rsb_file <- file.path(merged_dir, paste0(comp, ".genome_wide_rsb.txt.gz"))

  if (!file.exists(rsb_file)) {
    warning("  Rsb file not found: ", rsb_file, "\n")
    next
  }

  # Read genome-wide Rsb
  rsb_data <- read_tsv(rsb_file, show_col_types = FALSE)

  # Find Rsb column
  rsb_col <- grep("^RSB", names(rsb_data), value = TRUE)[1]

  # Extract values within each peak
  peak_rsb_list <- list()

  for (i in 1:nrow(peaks)) {
    peak_chr <- peaks$chromosome[i]
    peak_start <- peaks$start[i]
    peak_end <- peaks$end[i]
    peak_num <- peaks$peak[i]

    cat("  Extracting peak", peak_num, "(", peak_chr, ":", peak_start, "-", peak_end, ")... ")

    # Filter Rsb data to peak region
    peak_rsb <- rsb_data %>%
      filter(CHR == peak_chr,
             POSITION >= peak_start,
             POSITION <= peak_end) %>%
      mutate(peak = peak_num,
             comparison = comp,
             Rsb = .data[[rsb_col]])  # Rename to standard column

    cat(nrow(peak_rsb), "positions\n")

    peak_rsb_list[[i]] <- peak_rsb
  }

  # Combine all peaks for this comparison
  comp_peak_rsb <- bind_rows(peak_rsb_list)
  all_peak_rsb[[comp]] <- comp_peak_rsb

  cat("  Total peak positions:", nrow(comp_peak_rsb), "\n\n")
}

# Combine all comparisons
combined_peak_rsb <- bind_rows(all_peak_rsb)

cat("Total rows in combined Rsb data:", nrow(combined_peak_rsb), "\n")
cat("Comparisons:", paste(unique(combined_peak_rsb$comparison), collapse = ", "), "\n")
cat("Peaks:", paste(unique(combined_peak_rsb$peak), collapse = ", "), "\n\n")

# Save combined peak Rsb values
rsb_outfile <- file.path(outdir, "peak_rsb_values.txt")
write_tsv(combined_peak_rsb, rsb_outfile)
cat("Saved peak Rsb values:", rsb_outfile, "\n")
cat("  File size:", round(file.size(rsb_outfile) / 1024^2, 2), "MB\n\n")

# Calculate summary statistics per peak per comparison
peak_rsb_summary <- combined_peak_rsb %>%
  group_by(comparison, peak, CHR) %>%
  summarize(
    n_positions = n(),
    mean_Rsb = mean(Rsb, na.rm = TRUE),
    median_Rsb = median(Rsb, na.rm = TRUE),
    max_Rsb = max(Rsb, na.rm = TRUE),
    min_Rsb = min(Rsb, na.rm = TRUE),
    max_abs_Rsb = max(abs(Rsb), na.rm = TRUE),
    sd_Rsb = sd(Rsb, na.rm = TRUE),
    position_max_abs_Rsb = POSITION[which.max(abs(Rsb))],
    .groups = "drop"
  ) %>%
  arrange(desc(max_abs_Rsb))

# Save summary
rsb_summary_outfile <- file.path(outdir, "peak_rsb_summary.txt")
write_tsv(peak_rsb_summary, rsb_summary_outfile)
cat("Saved Rsb summary:", rsb_summary_outfile, "\n\n")

# ======================================================================
# Part 5: Calculate Percentile Rankings for Peak Maxima (inES/iHS)
# ======================================================================

cat("========================================\n")
cat("Part 5: Calculating Percentile Rankings (inES/iHS)\n")
cat("========================================\n\n")

ines_percentile_list <- list()

for (pop in pops) {
  cat("Processing population:", pop, "\n")

  ines_file <- file.path(merged_dir, paste0(pop, ".genome_wide_ines.txt.gz"))

  if (!file.exists(ines_file)) {
    warning("  inES file not found: ", ines_file, "\n")
    next
  }

  # Read genome-wide inES/iHS
  ines_data <- read_tsv(ines_file, show_col_types = FALSE)

  # Get peak summaries for this population
  pop_peak_summary <- peak_ines_summary %>%
    filter(population == pop)

  # Calculate percentiles for each peak
  for (i in 1:nrow(pop_peak_summary)) {
    peak_num <- pop_peak_summary$peak[i]
    chr <- pop_peak_summary$CHR[i]

    # inES percentile (max value)
    max_ines <- pop_peak_summary$max_inES[i]
    percentile_ines <- mean(ines_data$inES <= max_ines, na.rm = TRUE) * 100
    rank_ines <- sum(ines_data$inES > max_ines, na.rm = TRUE) + 1

    # iHS percentile (max absolute value)
    max_abs_ihs <- pop_peak_summary$max_abs_iHS[i]
    percentile_ihs <- mean(abs(ines_data$iHS) <= max_abs_ihs, na.rm = TRUE) * 100
    rank_ihs <- sum(abs(ines_data$iHS) > max_abs_ihs, na.rm = TRUE) + 1

    # Store results
    ines_percentile_list[[paste(pop, peak_num, "inES", sep = "_")]] <- tibble(
      population = pop,
      peak = peak_num,
      chromosome = chr,
      metric = "inES",
      peak_max_value = max_ines,
      percentile = percentile_ines,
      rank = rank_ines,
      total_snps = nrow(ines_data)
    )

    ines_percentile_list[[paste(pop, peak_num, "iHS", sep = "_")]] <- tibble(
      population = pop,
      peak = peak_num,
      chromosome = chr,
      metric = "iHS",
      peak_max_value = max_abs_ihs,
      percentile = percentile_ihs,
      rank = rank_ihs,
      total_snps = nrow(ines_data)
    )

    cat("  Peak", peak_num, "- inES:", round(percentile_ines, 2),
        "percentile | iHS:", round(percentile_ihs, 2), "percentile\n")
  }

  cat("\n")
}

# ======================================================================
# Part 6: Calculate Percentile Rankings for Peak Maxima (Rsb)
# ======================================================================

cat("========================================\n")
cat("Part 6: Calculating Percentile Rankings (Rsb)\n")
cat("========================================\n\n")

rsb_percentile_list <- list()

for (comp in comparisons) {
  cat("Processing comparison:", comp, "\n")

  rsb_file <- file.path(merged_dir, paste0(comp, ".genome_wide_rsb.txt.gz"))

  if (!file.exists(rsb_file)) {
    warning("  Rsb file not found: ", rsb_file, "\n")
    next
  }

  # Read genome-wide Rsb
  rsb_data <- read_tsv(rsb_file, show_col_types = FALSE)
  rsb_col <- grep("^RSB", names(rsb_data), value = TRUE)[1]

  # Get peak summaries for this comparison
  comp_peak_summary <- peak_rsb_summary %>%
    filter(comparison == comp)

  # Calculate percentiles for each peak
  for (i in 1:nrow(comp_peak_summary)) {
    peak_num <- comp_peak_summary$peak[i]
    chr <- comp_peak_summary$CHR[i]

    # Rsb percentile (max absolute value)
    max_abs_rsb <- comp_peak_summary$max_abs_Rsb[i]
    percentile_rsb <- mean(abs(rsb_data[[rsb_col]]) <= max_abs_rsb, na.rm = TRUE) * 100
    rank_rsb <- sum(abs(rsb_data[[rsb_col]]) > max_abs_rsb, na.rm = TRUE) + 1

    # Store results
    rsb_percentile_list[[paste(comp, peak_num, sep = "_")]] <- tibble(
      comparison = comp,
      peak = peak_num,
      chromosome = chr,
      metric = "Rsb",
      peak_max_value = max_abs_rsb,
      percentile = percentile_rsb,
      rank = rank_rsb,
      total_snps = nrow(rsb_data)
    )

    cat("  Peak", peak_num, "- Rsb:", round(percentile_rsb, 2), "percentile\n")
  }

  cat("\n")
}

# ======================================================================
# Part 7: Save Percentile Rankings
# ======================================================================

cat("========================================\n")
cat("Part 7: Saving Percentile Rankings\n")
cat("========================================\n\n")

# Combine inES/iHS percentiles
ines_percentiles <- bind_rows(ines_percentile_list) %>%
  arrange(desc(percentile))

# Combine Rsb percentiles
rsb_percentiles <- bind_rows(rsb_percentile_list) %>%
  arrange(desc(percentile))

# Save percentiles
ines_percentile_file <- file.path(outdir, "peak_ines_percentiles.txt")
write_tsv(ines_percentiles, ines_percentile_file)
cat("Saved inES/iHS percentiles:", ines_percentile_file, "\n")

rsb_percentile_file <- file.path(outdir, "peak_rsb_percentiles.txt")
write_tsv(rsb_percentiles, rsb_percentile_file)
cat("Saved Rsb percentiles:", rsb_percentile_file, "\n\n")

# Print top percentile rankings
cat("Top percentile rankings (inES/iHS):\n")
print(ines_percentiles %>% head(10), n = 10)

cat("\nTop percentile rankings (Rsb):\n")
print(rsb_percentiles %>% head(10), n = 10)

# ======================================================================
# Part 8: Calculate Empirical P-Values
# ======================================================================

cat("\n========================================\n")
cat("Part 8: Calculating Empirical P-Values\n")
cat("========================================\n\n")

ines_pvalue_list <- list()

for (pop in pops) {
  ines_file <- file.path(merged_dir, paste0(pop, ".genome_wide_ines.txt.gz"))

  if (!file.exists(ines_file)) {
    next
  }

  ines_data <- read_tsv(ines_file, show_col_types = FALSE)

  pop_peak_summary <- peak_ines_summary %>%
    filter(population == pop)

  for (i in 1:nrow(pop_peak_summary)) {
    peak_num <- pop_peak_summary$peak[i]
    chr <- pop_peak_summary$CHR[i]

    # inES empirical p-value (upper tail)
    max_ines <- pop_peak_summary$max_inES[i]
    p_upper_ines <- mean(ines_data$inES >= max_ines, na.rm = TRUE)
    p_two_tail_ines <- 2 * min(p_upper_ines, 1 - p_upper_ines)

    # iHS empirical p-value (two-tailed for absolute values)
    max_abs_ihs <- pop_peak_summary$max_abs_iHS[i]
    p_upper_ihs <- mean(abs(ines_data$iHS) >= max_abs_ihs, na.rm = TRUE)

    ines_pvalue_list[[paste(pop, peak_num, "inES", sep = "_")]] <- tibble(
      population = pop,
      peak = peak_num,
      chromosome = chr,
      metric = "inES",
      peak_max_value = max_ines,
      p_value_upper = p_upper_ines,
      p_value_two_tailed = p_two_tail_ines
    )

    ines_pvalue_list[[paste(pop, peak_num, "iHS", sep = "_")]] <- tibble(
      population = pop,
      peak = peak_num,
      chromosome = chr,
      metric = "iHS",
      peak_max_value = max_abs_ihs,
      p_value_upper = p_upper_ihs,
      p_value_two_tailed = NA_real_  # Two-tailed not meaningful for |iHS|
    )
  }
}

# Rsb p-values
rsb_pvalue_list <- list()

for (comp in comparisons) {
  rsb_file <- file.path(merged_dir, paste0(comp, ".genome_wide_rsb.txt.gz"))

  if (!file.exists(rsb_file)) {
    next
  }

  rsb_data <- read_tsv(rsb_file, show_col_types = FALSE)
  rsb_col <- grep("^RSB", names(rsb_data), value = TRUE)[1]

  comp_peak_summary <- peak_rsb_summary %>%
    filter(comparison == comp)

  for (i in 1:nrow(comp_peak_summary)) {
    peak_num <- comp_peak_summary$peak[i]
    chr <- comp_peak_summary$CHR[i]

    max_abs_rsb <- comp_peak_summary$max_abs_Rsb[i]
    p_upper_rsb <- mean(abs(rsb_data[[rsb_col]]) >= max_abs_rsb, na.rm = TRUE)

    rsb_pvalue_list[[paste(comp, peak_num, sep = "_")]] <- tibble(
      comparison = comp,
      peak = peak_num,
      chromosome = chr,
      metric = "Rsb",
      peak_max_value = max_abs_rsb,
      p_value_upper = p_upper_rsb,
      p_value_two_tailed = NA_real_
    )
  }
}

# Combine and save p-values
ines_pvalues <- bind_rows(ines_pvalue_list) %>%
  arrange(p_value_upper)

rsb_pvalues <- bind_rows(rsb_pvalue_list) %>%
  arrange(p_value_upper)

ines_pvalue_file <- file.path(outdir, "peak_ines_empirical_pvalues.txt")
write_tsv(ines_pvalues, ines_pvalue_file)
cat("Saved inES/iHS p-values:", ines_pvalue_file, "\n")

rsb_pvalue_file <- file.path(outdir, "peak_rsb_empirical_pvalues.txt")
write_tsv(rsb_pvalues, rsb_pvalue_file)
cat("Saved Rsb p-values:", rsb_pvalue_file, "\n\n")

# ======================================================================
# Summary
# ======================================================================

cat("\n========================================\n")
cat("SUMMARY\n")
cat("========================================\n")
cat("Genome-wide summary statistics:\n")
cat("  inES/iHS:", ines_summary_file, "\n")
cat("  Rsb:", rsb_summary_file, "\n\n")

cat("Peak extractions:\n")
cat("  inES/iHS values:", ines_outfile, "(", nrow(combined_peak_ines), "positions )\n")
cat("  Rsb values:", rsb_outfile, "(", nrow(combined_peak_rsb), "positions )\n")
cat("  inES/iHS summary:", ines_summary_outfile, "\n")
cat("  Rsb summary:", rsb_summary_outfile, "\n\n")

cat("Percentile rankings:\n")
cat("  inES/iHS:", ines_percentile_file, "\n")
cat("  Rsb:", rsb_percentile_file, "\n\n")

cat("Empirical p-values:\n")
cat("  inES/iHS:", ines_pvalue_file, "\n")
cat("  Rsb:", rsb_pvalue_file, "\n\n")

cat("========================================\n")
cat("Peak EHH extraction complete\n")
cat("========================================\n")
cat("\nNext step: bash code/12_Sweep_Detection/ehh/05_render_ehh_report.sh\n")
