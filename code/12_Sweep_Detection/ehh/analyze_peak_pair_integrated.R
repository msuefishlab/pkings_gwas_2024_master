library(ggplot2)
library(patchwork)

# helpers --------------------------------------------------------------

.find_col <- function(nms, candidates) {
  hit <- dplyr::first(intersect(candidates, nms))
  if (is.na(hit)) stop("Could not find any of: ", paste(candidates, collapse=", "))
  hit
}

library(readr)
library(dplyr)

# Strip ".1/.2" or "_1/_2" from hap names to get the individual ID
.ind_from_hap <- function(haps) sub("(\\.|_)?[12]$", "", haps)

# Read ID lists from files (or pass character vectors directly)
read_id_list <- function(path_or_vec) {
  if (length(path_or_vec) == 1 && file.exists(path_or_vec)) {
    readr::read_table(path_or_vec, col_names = FALSE, show_col_types = FALSE)[[1]]
  } else {
    as.character(path_or_vec)
  }
}

# ---- MAIN SPLITTER -------------------------------------------------------
# Returns: list with (1) subsets hh_BP / hh_TP, (2) indices, and (3) a named tibble
split_hh_by_groups <- function(hh_all, bp_ids, tp_ids, keep_unassigned = FALSE) {
  # normalize inputs
  bp_ids <- read_id_list(bp_ids)
  tp_ids <- read_id_list(tp_ids)
  
  haps <- hap.names(hh_all)
  indiv <- .ind_from_hap(haps)
  
  # map each haplotype to BP/TP/other
  group <- case_when(
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
  map_df <- tibble(
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


# core analysis --------------------------------------------------------

#' Analyze one peak for a BP/TP pair: EHHS@core, inES (per pop), Rsb (inES-based)
#' @return list(meta, df = list(ehhs, ines, rsb))
analyze_peak_tang <- function(hh_bp, hh_tp, peak_info, labels = c("BP","TP"), win = 5e5) {
  stopifnot(nrow(peak_info) == 1, all(c("chrom","idx_bp") %in% names(peak_info)))
  
  # nearest markers to the provided genomic coordinate
  idx_bp <- nearest_idx(hh_bp, peak_info$idx_bp)
  idx_tp <- nearest_idx(hh_tp, peak_info$idx_bp)
  pos_bp <- get_pos(hh_bp)[idx_bp]
  pos_tp <- get_pos(hh_tp)[idx_tp]
  
  rng <- c(max(0, peak_info$idx_bp - win), peak_info$idx_bp + win)
  
  # EHHS at the core SNP (both pops)
  ehhs_bp <- calc_ehhs(hh_bp, mrk = idx_bp, include_zero_values = TRUE)
  ehhs_tp <- calc_ehhs(hh_tp, mrk = idx_tp, include_zero_values = TRUE)
  
  # Tidy EHHS (be forgiving about column names)
  as_ehhs_df <- function(obj, poplab) {
    df <- tibble::as_tibble(obj$ehhs %||% obj$EHH)
    poscol <- .find_col(names(df), c("POSITION","position","POS","pos"))
    valcols <- grep("EHH", names(df), value = TRUE)
    if (length(valcols) == 0) {
      valcols <- setdiff(names(df)[sapply(df, is.numeric)], poscol)
    }
    df |>
      dplyr::select(dplyr::all_of(c(poscol, valcols))) |>
      dplyr::rename(POSITION = !!poscol) |>
      tidyr::pivot_longer(dplyr::all_of(valcols), names_to = "metric", values_to = "value") |>
      dplyr::mutate(POP = poplab)
  }
  
  df_ehhs <- dplyr::bind_rows(
    as_ehhs_df(ehhs_bp, labels[1]),
    as_ehhs_df(ehhs_tp, labels[2])
  ) |>
    dplyr::filter(dplyr::between(POSITION, rng[1], rng[2]))
  
  # inES per pop (via scan_hh)
  scan_bp <- scan_hh(hh_bp, polarized = TRUE, discard_integration_at_border = TRUE)
  scan_tp <- scan_hh(hh_tp, polarized = TRUE, discard_integration_at_border = TRUE)
  
  to_ines_df <- function(tab, poplab, chrom, rng) {
    df <- tibble::as_tibble(tab$scan_hh %||% tab)  # support older/newer rehh
    chrcol <- .find_col(names(df), c("CHR","chr","CHROM","chrom"))
    poscol <- .find_col(names(df), c("POSITION","position","POS","pos"))
    incol  <- .find_col(names(df), c("inES","INES","ines"))
    df |>
      dplyr::filter(.data[[chrcol]] == chrom,
                    dplyr::between(.data[[poscol]], rng[1], rng[2])) |>
      dplyr::transmute(CHR = .data[[chrcol]],
                       POSITION = .data[[poscol]],
                       inES = .data[[incol]],
                       POP = poplab)
  }
  
  df_ines <- dplyr::bind_rows(
    to_ines_df(scan_bp, labels[1], peak_info$chrom, rng),
    to_ines_df(scan_tp, labels[2], peak_info$chrom, rng)
  )
  
  # Rsb from inES (Tang 2007): use ines2rsb
  rsb_tab <- ines2rsb(scan_pop1 = scan_bp, scan_pop2 = scan_tp,
                      popname1 = labels[1], popname2 = labels[2])
  expected <- paste0("RSB_", labels[1], "_", labels[2])
  rsbcol <- if (expected %in% names(rsb_tab)) {
    expected
  } else {
    # fall back: first column starting with "RSB"
    dplyr::first(grep("^RSB", names(rsb_tab), value = TRUE))
  }
  if (is.na(rsbcol)) stop("Could not find RSB column in rsb_tab")
  chrcol <- .find_col(names(rsb_tab), c("CHR","chr","CHROM","chrom"))
  poscol <- .find_col(names(rsb_tab), c("POSITION","position","POS","pos"))
  
  df_rsb <- tibble::as_tibble(rsb_tab) |>
    dplyr::filter(.data[[chrcol]] == peak_info$chrom,
                  dplyr::between(.data[[poscol]], rng[1], rng[2])) |>
    dplyr::transmute(CHR = .data[[chrcol]],
                     POSITION = .data[[poscol]],
                     Rsb = .data[[rsbcol]])
  
  list(
    meta = list(
      chrom   = peak_info$chrom,
      idx_bp  = peak_info$idx_bp,
      pos_bp  = pos_bp,
      pos_tp  = pos_tp,
      labels  = labels,
      window  = rng
    ),
    df = list(
      ehhs = df_ehhs,   # POSITION, metric, value, POP
      ines = df_ines,   # CHR, POSITION, inES, POP
      rsb  = df_rsb     # CHR, POSITION, Rsb
    )
  )
}




library(ggplot2)
library(patchwork)

#' Plot Tang-style 3-panel summary for a single peak result
#' Panel A: EHHS (site-specific, not normalized)
#' Panel B: inES (per pop)
#' Panel C: Rsb (between pops)
#' BP phenotypes always = #1ea8e0, TP phenotypes always = #ed2224
plot_peak_tang <- function(res) {
  chrom  <- res$meta$chrom
  center <- res$meta$idx_bp
  rng    <- res$meta$window
  labs2  <- res$meta$labels
  
  # colors by convention
  colmap <- c("BP" = "#1ea8e0", "TP" = "#ed2224")
  pop_colors <- setNames(
    ifelse(grepl("^BP", labs2), colmap["BP"], colmap["TP"]),
    labs2
  )
  
  # Centered coordinates
  ehhs_df <- res$df$ehhs |>
    dplyr::mutate(X = POSITION - center) |>
    dplyr::filter(metric %in% c("EHHS","ehhs","EHHS_pop","EHHS_All"))  # keep only EHHS-type
  
  ines_df <- res$df$ines |>
    dplyr::mutate(X = POSITION - center)
  
  rsb_df  <- res$df$rsb |>
    dplyr::mutate(X = POSITION - center)
  
  # Panel A: EHHS
  pA <- ggplot(ehhs_df, aes(X, value, color = POP)) +
    geom_line(na.rm = TRUE, linewidth = 0.8) +
    scale_color_manual(values = pop_colors) +
    geom_vline(xintercept = 0, linetype = 3) +
    coord_cartesian(xlim = c(rng[1] - center, rng[2] - center)) +
    labs(title = "A. EHHS at core SNP",
         x = "Distance from core SNP (bp)", y = "EHHS") +
    theme_bw()
  
  # Panel B: inES
  pB <- ggplot(ines_df, aes(X, inES, color = POP)) +
    geom_point(size = 0.8, alpha = 0.9, na.rm = TRUE) +
    geom_smooth(se = FALSE, method = "loess", formula = y ~ x, span = 0.3) +
    scale_color_manual(values = pop_colors) +
    geom_vline(xintercept = 0, linetype = 3) +
    coord_cartesian(xlim = c(rng[1] - center, rng[2] - center)) +
    labs(title = sprintf("B. inES (%s vs %s)", labs2[1], labs2[2]),
         x = "Distance from core SNP (bp)", y = "inES") +
    theme_bw()
  
  # Panel C: Rsb
  pC <- ggplot(rsb_df, aes(X, Rsb)) +
    geom_hline(yintercept = c(-2, 2), linetype = 2) +
    geom_point(size = 0.8, alpha = 0.9, na.rm = TRUE, color = "black") +
    geom_vline(xintercept = 0, linetype = 3) +
    coord_cartesian(xlim = c(rng[1] - center, rng[2] - center)) +
    labs(title = sprintf("C. Rsb (%s/%s)", labs2[1], labs2[2]),
         x = "Distance from core SNP (bp)", y = "Rsb") +
    theme_bw()
  
  (pA / pB / pC) +
    plot_annotation(
      title = sprintf("Peak near %s:%s — %s vs %s",
                      chrom, center, labs2[1], labs2[2])
    )
}
