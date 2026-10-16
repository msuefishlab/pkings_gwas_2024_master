#!/usr/bin/env Rscript
## 12_tail_figure.R: the tail figure of DESIGN.md section 10. For each pair and analysis, the genome-wide count
## of sites at least as extreme as delta, observed outside every GWAS peak and expected under each admissible
## model (tail/curve.tsv from 09_tail.py), with each pair's own index SNPs marked at their tail threshold in
## the analysis that is their headline (tail/unit_results.tsv).
##   Rscript code/11_Demographic_Modeling/fixed_differences/12_tail_figure.R <tail_dir> <out_prefix>
## writes <out_prefix>.pdf and <out_prefix>.png
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
tail_dir <- args[1]
out <- args[2]

FLOOR <- 0.01   # bottom of the log axis: observed zeros are marked there, model curves stop there
pair_lab <- c(A1w = "BAM vs APA + BEN", A2w = "BIK vs DOV + DOG", A3w = "BAVA vs DOV + DOG",
              A4w = "BIR vs DOV + DOG", C1w = "APA vs BEN (control)", C2w = "DOV vs DOG (control)")
short_lab <- c(A1w = "BAM", A2w = "BIK", A3w = "BAVA", A4w = "BIR", C1w = "APA-BEN", C2w = "DOV-DOG")
an_lab <- c(unfolded = "Polarized (BP-derived allele)", folded = "Label-free (either allele)")

curve <- read_tsv(file.path(tail_dir, "curve.tsv"), show_col_types = FALSE)
units <- read_tsv(file.path(tail_dir, "unit_results.tsv"), show_col_types = FALSE,
                  col_types = cols(.default = col_character()))

best <- units |>
  filter(best_model != "none") |>
  distinct(pair, analysis, best_model)

d <- curve |>
  mutate(model = sub("f$", "", series),
         kind = ifelse(series == "observed_nonpeak", "observed", "model")) |>
  left_join(best, by = c("pair", "analysis")) |>
  mutate(is_best = kind == "model" & series == best_model,
         y = pmax(count_genome, FLOOR),
         pair = factor(pair_lab[pair], levels = pair_lab),
         analysis = factor(an_lab[analysis], levels = an_lab))

own <- units |>
  filter(own == "True", headline == "True") |>
  mutate(x = abs(as.numeric(tail_threshold)),
         O = pmax(as.numeric(O_genome), FLOOR),
         E = as.numeric(E_genome_best),
         model = sub("f$", "", best_model),
         pair_id = pair,
         pair = factor(pair_lab[pair], levels = pair_lab),
         analysis = factor(an_lab[analysis], levels = an_lab)) |>
  filter(x >= 0.3)

# plotted heights per row, so marker labels sit just above each row's data (the y range is free per row)
row_top <- d |>
  filter(kind == "observed" | count_genome >= FLOOR) |>
  group_by(pair) |>
  summarise(top = max(y), .groups = "drop")

marks <- own |>
  mutate(label = paste0("P", peak)) |>
  group_by(pair, analysis, x) |>
  summarise(label = paste(label, collapse = ","), .groups = "drop") |>
  left_join(row_top, by = "pair") |>
  # stagger nearby markers: alternate label height, and put every second label right of its line unless the
  # line is at the right edge
  arrange(pair, analysis, x) |>
  group_by(pair, analysis) |>
  mutate(k = seq_along(x),
         ly = ifelse(k %% 2 == 1, top * 30, top * 4),
         hj = ifelse(k %% 2 == 0 & x < 0.95, -0.08, 1.08)) |>
  ungroup()

# Okabe-Ito colours and one linetype per model, so models stay apart in greyscale and for colour-blind readers
model_cols <- c(M0 = "#009E73", M1 = "#E69F00", M2 = "#0072B2", M3 = "#CC79A7")
model_lty <- c(M0 = "dotted", M1 = "dashed", M2 = "solid", M3 = "longdash")
y_lab <- function(b) c("0.01", "0.1", "1", "10", "100", "1k", "10k", "100k", "1M", "10M", "100M",
                       "1G")[round(log10(b)) + 3]
# panels with no admissible model, named in the caption
none_adm <- units |>
  group_by(pair, analysis) |>
  summarise(none = all(best_model == "none"), .groups = "drop") |>
  filter(none) |>
  mutate(lab = paste0(short_lab[pair], ifelse(analysis == "folded", " label-free", " polarized")))
none_txt <- if (nrow(none_adm) == 0) "every panel has at least one" else
  paste0("none admissible for ", paste(none_adm$lab, collapse = ", "))
obs <- filter(d, kind == "observed", count_genome > 0)
# where the observed count first reaches zero: one open circle on the floor
obs0 <- d |>
  filter(kind == "observed", count_genome == 0) |>
  group_by(pair, analysis) |>
  slice_min(delta, n = 1) |>
  ungroup()
# model curves end where they fall below the floor rather than running along it
mod <- filter(d, kind == "model", count_genome >= FLOOR)

p <- ggplot() +
  geom_vline(data = marks, aes(xintercept = x), colour = "grey55", linetype = "dashed", linewidth = 0.35) +
  geom_text(data = marks, aes(x = x, y = ly, label = label, hjust = hj), size = 2.6, vjust = 1,
            colour = "grey25") +
  geom_line(data = filter(mod, !is_best), aes(delta, y, colour = model, linetype = model, group = series),
            linewidth = 0.45) +
  geom_line(data = filter(mod, is_best), aes(delta, y, colour = model, linetype = model, group = series),
            linewidth = 1.2) +
  geom_line(data = obs, aes(delta, y), linewidth = 0.5, colour = "black") +
  geom_point(data = obs0, aes(delta, y), shape = 1, size = 1.8, stroke = 0.6, colour = "black") +
  # each own unit's O (black circle) and E under the best model (diamond), as in unit_results.tsv
  geom_point(data = filter(own, model != "none", E >= FLOOR), aes(x, E, fill = model), shape = 23, size = 2,
             stroke = 0.4, colour = "black", show.legend = FALSE) +
  geom_point(data = own, aes(x, O), shape = 21, size = 1.6, fill = "black", colour = "white", stroke = 0.3) +
  scale_colour_manual(values = model_cols, name = "Model") +
  scale_fill_manual(values = model_cols, guide = "none") +
  scale_linetype_manual(values = model_lty, name = "Model") +
  scale_y_log10(breaks = 10^seq(-2, 9, 1), labels = y_lab) +
  scale_x_continuous(limits = c(0.3, 1), breaks = seq(0.3, 1, 0.1), expand = expansion(mult = c(0.01, 0.04))) +
  facet_grid(pair ~ analysis, scales = "free_y", labeller = labeller(pair = label_wrap_gen(14))) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 0.8)),
         linetype = guide_legend(nrow = 1)) +
  labs(x = expression("Allele-frequency difference at least " * Delta),
       y = "Sites in the genome (log scale)",
       caption = paste(strwrap(paste0(
         "Black line: observed outside the GWAS peaks; open circle: first threshold with no such site. Coloured ",
         "lines: expected under each admissible model, thick = best model (", none_txt, "). Dashed: own index ",
         "SNPs (P1-P6) at their tail threshold, in their headline analysis; at each, black dot = observed count ",
         "O and diamond = expected count E under the best model, as in the results table. Observed zeros are ",
         "drawn at 0.01."), 125), collapse = "\n")) +
  theme_bw(base_size = 9) +
  theme(strip.text.y = element_text(angle = 0, hjust = 0), panel.grid.minor = element_blank(),
        legend.position = "bottom", legend.title = element_text(size = 7.5),
        plot.caption = element_text(hjust = 0, size = 7.5))

ggsave(paste0(out, ".pdf"), p, width = 7.2, height = 10.5)
ggsave(paste0(out, ".png"), p, width = 7.2, height = 10.5, dpi = 200)
cat("wrote", paste0(out, c(".pdf", ".png")), "\n")
