#!/usr/bin/env Rscript
## 09_fig_s7_revised.R
## Rebuild Fig. S7 (EHHS / inES / Rsb at the six GWAS peaks) with Peak 6 shown for BP2 vs TP2
##
## The published Fig. S7 was assembled from the final chunk of ehh_merge_pheno.Rmd
## (objects p1.1, p2.2, p3.1, p4.1, p5.1, p6.2 -> peak_ehh_stats.svg, then edited by hand).
## Peak 6 used BP3 vs TP2 (p6.2), but the Peak 6 BP haplotype is BP2-specific. This script
## redraws the same six columns with p6.1 (BP2 vs TP2) in place of p6.2. Everything else is
## unchanged: the same haplotype objects, the same per-peak computation
## (analyze_peak_tang() in analyze_peak_pair_integrated.R), +/-250 kb around the index SNP,
## the same y-limits and the same genome-wide 95th/99th percentile reference lines.
##
## What the panels show (unchanged from the published figure):
##   Row 1  EHHS: site-specific EHH (rehh::calc_ehhs) at the index SNP, one curve per
##          POPULATION (BP vs TP haplotypes). Not allele-specific EHH.
##   Row 2  inES per SNP (rehh::scan_hh on the whole chromosome, per population); dashed lines
##          are that population's genome-wide 95th (thin) and 99th (thick) inES percentiles.
##   Row 3  Rsb (rehh::ines2rsb, BP vs TP, standardized within the chromosome); grey lines at
##          +/-2, orange/red dashed lines at +/- the genome-wide 95th/99th Rsb percentiles.
##
## Usage:
##   Rscript code/12_Sweep_Detection/ehh/09_fig_s7_revised.R
##
## Prerequisites (local, no HPCC needed):
##   - output_data/12_Sweep_Detection/processed_rehh_data.Rdata (haplohh objects; built by the
##     commented-out chunk of ehh_merge_pheno.Rmd from output_data/12_Sweep_Detection/by_chrom/)
##   - output_data/12_Sweep_Detection/ehh/genome_wide_{ines,rsb}_summary.txt (06_extract_peak_ehh.R)
##   - input_data/10_Introgression/peaks.bed (extended peaks; index SNPs are identical to the
##     core PEAKS.bed used by the Rmd, see decisions D-002)
##
## Outputs: output_data/12_Sweep_Detection/ehh/fig_s7_revised/
##   fig_s7_revised.pdf / .png         six columns, Peak 6 = BP2 vs TP2
##   fig_s7_peak6_bp3_vs_bp2.pdf / .png  Peak 6 for BP3 (published) and BP2 side by side
##   fig_s7_peak_results.rds           cached analyze_peak_tang() results (delete to recompute)
##
## Author: Jason Gallant Lab
## Date: 2026

suppressPackageStartupMessages({
  library(tidyverse)
  library(rprojroot)
  library(rehh)
  library(patchwork)
})

root <- find_root(".git/index")

cat("========================================\n")
cat("Fig. S7 (revised): Peak 6 as BP2 vs TP2\n")
cat("========================================\n\n")

outdir <- file.path(root, "output_data/12_Sweep_Detection/ehh/fig_s7_revised")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Helpers that ehh_merge_pheno.Rmd defines in its setup chunk (needed by analyze_peak_tang)
get_pos <- function(hh) {
  if (!is.null(hh@positions)) return(hh@positions)
  if (!is.null(hh@position))  return(hh@position)
  stop("No positions/position slot found on haplohh.")
}
nearest_idx <- function(hh, pos) {
  pos_vec <- as.numeric(get_pos(hh))
  which.min(abs(pos_vec - pos))
}

# analyze_peak_tang() is unchanged from the committed version (only plot_peak_tang differs)
source(file.path(root, "code/12_Sweep_Detection/ehh/analyze_peak_pair_integrated.R"))

# ----------------------------------------------------------------------
# Inputs
# ----------------------------------------------------------------------
peaks <- read.table(file.path(root, "input_data/10_Introgression/peaks.bed"),
                    col.names = c("chrom", "start", "end", "peak", "idx_p", "idx_bp"))
peaks$chrom <- paste0("chr", peaks$chrom)

ines_summary <- read_tsv(file.path(root, "output_data/12_Sweep_Detection/ehh/genome_wide_ines_summary.txt"),
                         show_col_types = FALSE)
rsb_summary  <- read_tsv(file.path(root, "output_data/12_Sweep_Detection/ehh/genome_wide_rsb_summary.txt"),
                         show_col_types = FALSE)

# Panels: published columns, except Peak 6 now BP2 (published used BP3).
# peak6_BP3 is computed too, for the side-by-side comparison figure.
panels <- tribble(
  ~key,        ~peak, ~obj,              ~labels,
  "peak1_BP1", 1,     "peak1.BP1.splt",  c("BP1", "TP1"),
  "peak2_BP3", 2,     "peak2.BP3.splt",  c("BP3", "TP2"),
  "peak3_BP2", 3,     "peak3.BP2.splt",  c("BP2", "TP2"),
  "peak4_BP1", 4,     "peak4.BP1.splt",  c("BP1", "TP1"),
  "peak5_BP1", 5,     "peak5.BP1.splt",  c("BP1", "TP1"),
  "peak6_BP2", 6,     "peak6.BP2.splt",  c("BP2", "TP2"),
  "peak6_BP3", 6,     "peak6.BP3.splt",  c("BP3", "TP2")
)

cache_file <- file.path(outdir, "fig_s7_peak_results.rds")
if (file.exists(cache_file)) {
  cat("Loading cached peak results:", cache_file, "\n\n")
  res <- readRDS(cache_file)
} else {
  rdata <- file.path(root, "output_data/12_Sweep_Detection/processed_rehh_data.Rdata")
  if (!file.exists(rdata)) stop("Missing haplotype objects: ", rdata)
  cat("Loading haplotype objects:", rdata, "\n")
  hh_env <- new.env()
  load(rdata, envir = hh_env)

  res <- list()
  for (i in seq_len(nrow(panels))) {
    p <- panels[i, ]
    cat("  ", p$key, ": scan_hh / calc_ehhs / ines2rsb on", peaks$chrom[peaks$peak == p$peak], "... ")
    splt <- get(p$obj, envir = hh_env)
    t0 <- Sys.time()
    res[[p$key]] <- analyze_peak_tang(splt$hh$BP, splt$hh$TP, peaks[peaks$peak == p$peak, ],
                                      labels = p$labels[[1]], win = 2.5e5)
    cat(round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
  }
  saveRDS(res, cache_file)
  cat("Saved cache:", cache_file, "\n\n")
}

# ----------------------------------------------------------------------
# Plotting (layout of the published figure: kb axis, legends inside panels)
# ----------------------------------------------------------------------
colmap <- c(BP = "#1ea8e0", TP = "#ed2224")
ines_ylim <- c(0, 60000)
rsb_ylim  <- c(-3, 4)
# Sized for a full-width ESM page (published image is ~1.1:1); fonts are legible at this size
FIG_W <- 7.5; FIG_H <- 6.7; BASE <- 8; PT <- 0.25

plot_column <- function(r, title, show_y = TRUE) {
  labs2  <- r$meta$labels
  center <- r$meta$idx_bp
  xlim   <- (r$meta$window - center) / 1000
  pop_cols <- setNames(ifelse(grepl("^BP", labs2), colmap["BP"], colmap["TP"]), labs2)

  # TP is drawn first so the BP curve sits on top; the legend still lists BP first
  ehhs_df <- r$df$ehhs |>
    filter(metric %in% c("EHHS", "ehhs", "EHHS_pop", "EHHS_All")) |>
    mutate(X = (POSITION - center) / 1000, POP = factor(POP, levels = rev(labs2))) |>
    arrange(POP)
  ines_df <- r$df$ines |> mutate(X = (POSITION - center) / 1000, POP = factor(POP, levels = rev(labs2))) |>
    arrange(POP)
  rsb_df  <- r$df$rsb  |> mutate(X = (POSITION - center) / 1000)

  base_theme <- theme_bw(base_size = BASE) +
    theme(legend.position = "inside",
          legend.position.inside = c(0, 1),
          legend.margin = margin(1, 1, 1, 1),
          legend.justification = c(0, 1),
          legend.title = element_blank(),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.key.size = unit(0.6, "lines"),
          legend.key.width = unit(0.5, "lines"),
          legend.text = element_text(size = BASE - 1.5, margin = margin(l = 1.5)),
          axis.title.x = element_blank(),
          axis.text.x = element_blank())

  pA <- ggplot(ehhs_df, aes(X, value, color = POP)) +
    geom_line(na.rm = TRUE, linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = 3) +
    scale_color_manual(values = pop_cols, breaks = labs2) +
    coord_cartesian(xlim = xlim, ylim = c(0, 1)) +
    labs(title = title, y = "EHHS") +
    base_theme +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = BASE + 2))

  pB <- ggplot(ines_df, aes(X, inES, color = POP)) +
    geom_point(size = PT, alpha = 0.9, na.rm = TRUE, show.legend = FALSE) +
    geom_smooth(se = FALSE, method = "loess", formula = y ~ x, span = 0.3, na.rm = TRUE, linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = 3) +
    scale_color_manual(values = pop_cols, breaks = labs2) +
    coord_cartesian(xlim = xlim, ylim = ines_ylim) +
    labs(y = "inES") +
    base_theme
  for (pop in labs2) {
    s <- ines_summary |> filter(population == pop, metric == "inES")
    pB <- pB +
      geom_hline(yintercept = s$p95, linetype = "dashed", color = pop_cols[pop], alpha = 0.7, linewidth = 0.4) +
      geom_hline(yintercept = s$p99, linetype = "dashed", color = pop_cols[pop], alpha = 0.8, linewidth = 0.45)
  }

  rs <- rsb_summary |> filter(comparison == paste0(labs2[1], "_", labs2[2]))
  pC <- ggplot(rsb_df, aes(X, Rsb)) +
    geom_hline(yintercept = c(-2, 2), linetype = 2, color = "gray40", alpha = 0.6, linewidth = 0.4) +
    geom_hline(yintercept = c(-1, 1) * rs$p95, linetype = "dashed", color = "darkorange", alpha = 0.8, linewidth = 0.4) +
    geom_hline(yintercept = c(-1, 1) * rs$p99, linetype = "dashed", color = "red", alpha = 0.7, linewidth = 0.45) +
    geom_point(size = PT, alpha = 0.9, na.rm = TRUE, color = "black") +
    geom_vline(xintercept = 0, linetype = 3) +
    coord_cartesian(xlim = xlim, ylim = rsb_ylim) +
    scale_x_continuous(breaks = c(-200, -100, 0, 100, 200)) +
    labs(y = "Rsb", x = "Distance from index SNP (kb)") +
    theme_bw(base_size = BASE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = BASE - 1.5),
          axis.title.x = element_blank())

  if (!show_y) {
    pA <- pA + theme(axis.title.y = element_blank())
    pB <- pB + theme(axis.title.y = element_blank())
    pC <- pC + theme(axis.title.y = element_blank())
  }
  pA / pB / pC
}

assemble <- function(cols) {
  wrap_plots(cols, nrow = 1) +
    plot_annotation(caption = "Distance from index SNP (kb)",
                    theme = theme(plot.caption = element_text(hjust = 0.5, size = BASE + 1),
                                  plot.background = element_rect(fill = "white", colour = NA)))
}

order_keys <- c("peak1_BP1", "peak2_BP3", "peak3_BP2", "peak4_BP1", "peak5_BP1", "peak6_BP2")
cols <- map2(order_keys, seq_along(order_keys),
             ~ plot_column(res[[.x]], paste("Peak", panels$peak[panels$key == .x]), show_y = .y == 1))
fig <- assemble(cols)

ggsave(file.path(outdir, "fig_s7_revised.pdf"), fig, width = FIG_W, height = FIG_H, units = "in")
ggsave(file.path(outdir, "fig_s7_revised.png"), fig, width = FIG_W, height = FIG_H, units = "in", dpi = 600, bg = "white")
cat("Saved:", file.path(outdir, "fig_s7_revised.{pdf,png}"), "\n")

# Peak 6: published BP3 column next to the corrected BP2 column
cmp <- assemble(list(plot_column(res[["peak6_BP3"]], "Peak 6: BP3 (old)", show_y = TRUE),
                     plot_column(res[["peak6_BP2"]], "Peak 6: BP2 (new)", show_y = FALSE)))
ggsave(file.path(outdir, "fig_s7_peak6_bp3_vs_bp2.pdf"), cmp, width = 3.4, height = FIG_H, units = "in")
ggsave(file.path(outdir, "fig_s7_peak6_bp3_vs_bp2.png"), cmp, width = 3.4, height = FIG_H, units = "in", dpi = 600, bg = "white")
cat("Saved:", file.path(outdir, "fig_s7_peak6_bp3_vs_bp2.{pdf,png}"), "\n\n")

# Per-panel values at the index SNP, for the legend
idx_tab <- imap_dfr(res, function(r, key) {
  ip <- r$meta$idx_bp
  ines_at <- r$df$ines |> filter(POSITION == ip)
  tibble(panel = key, chrom = r$meta$chrom, index_snp = ip,
         nearest_marker_BP = r$meta$pos_bp, nearest_marker_TP = r$meta$pos_tp,
         inES_BP_at_index = ines_at$inES[ines_at$POP == r$meta$labels[1]][1],
         inES_TP_at_index = ines_at$inES[ines_at$POP == r$meta$labels[2]][1],
         Rsb_at_index = r$df$rsb$Rsb[r$df$rsb$POSITION == ip][1],
         max_Rsb_pm250kb = max(r$df$rsb$Rsb, na.rm = TRUE),
         min_Rsb_pm250kb = min(r$df$rsb$Rsb, na.rm = TRUE),
         pos_max_Rsb_pm250kb = r$df$rsb$POSITION[which.max(r$df$rsb$Rsb)])
})
write_tsv(idx_tab, file.path(outdir, "fig_s7_index_snp_values.tsv"))
print(idx_tab, width = Inf)

cat("\n========================================\n")
cat("Fig. S7 revision complete\n")
cat("========================================\n")
