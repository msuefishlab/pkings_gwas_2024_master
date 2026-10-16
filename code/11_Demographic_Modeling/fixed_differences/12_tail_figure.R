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

marks <- units |>
  filter(own == "True", headline == "True") |>
  mutate(x = abs(as.numeric(tail_threshold)),
         pair = factor(pair_lab[pair], levels = pair_lab),
         analysis = factor(an_lab[analysis], levels = an_lab),
         label = paste0("P", peak)) |>
  filter(x >= 0.3) |>
  group_by(pair, analysis, x) |>
  summarise(label = paste(label, collapse = ","), .groups = "drop") |>
  # stagger labels of nearby markers
  arrange(pair, analysis, x) |>
  group_by(pair, analysis) |>
  mutate(ly = ifelse(seq_along(x) %% 2 == 1, 3e7, 2e6)) |>
  ungroup()

model_cols <- c(M0 = "#1b9e77", M1 = "#d95f02", M2 = "#7570b3", M3 = "#e7298a")
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
  geom_text(data = marks, aes(x = x, y = ly, label = label), size = 2.6, hjust = 1.08, vjust = 1,
            colour = "grey25") +
  geom_line(data = filter(mod, !is_best), aes(delta, y, colour = model, group = series),
            linewidth = 0.4, alpha = 0.8) +
  geom_line(data = filter(mod, is_best), aes(delta, y, colour = model, group = series), linewidth = 1.1) +
  geom_line(data = obs, aes(delta, y), linewidth = 0.6, colour = "black") +
  geom_point(data = obs0, aes(delta, y), shape = 1, size = 1.8, stroke = 0.6, colour = "black") +
  scale_colour_manual(values = model_cols, name = "Model (admissible;\nthick = best)") +
  scale_y_log10(limits = c(FLOOR, 3e7), breaks = 10^seq(-2, 7, 1),
                labels = c("0.01", "0.1", "1", "10", "100", "1k", "10k", "100k", "1M", "10M")) +
  scale_x_continuous(limits = c(0.3, 1), breaks = seq(0.3, 1, 0.1)) +
  facet_grid(pair ~ analysis, labeller = labeller(pair = label_wrap_gen(10))) +
  labs(x = expression("Allele-frequency difference at least " * Delta),
       y = "Sites in the genome (log scale)",
       caption = paste(strwrap(paste0(
         "Black: observed outside the GWAS peaks; open circle: first threshold with no such site. Colour: ",
         "expected under each admissible model (", none_txt, "). Dashed: own index SNPs (P1-P6) at their tail ",
         "threshold, in their headline analysis."), 125), collapse = "\n")) +
  theme_bw(base_size = 9) +
  theme(strip.text.y = element_text(angle = 0, hjust = 0), panel.grid.minor = element_blank(),
        legend.position = "bottom", plot.caption = element_text(hjust = 0, size = 7.5))

ggsave(paste0(out, ".pdf"), p, width = 7.2, height = 10.5)
ggsave(paste0(out, ".png"), p, width = 7.2, height = 10.5, dpi = 200)
cat("wrote", paste0(out, c(".pdf", ".png")), "\n")
