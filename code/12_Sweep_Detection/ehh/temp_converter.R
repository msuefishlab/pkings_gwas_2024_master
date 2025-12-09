get_center_pos <- function(hh_obj, mrk_index) {
  # hh_obj uses genomic coordinates in @positions
  hh_obj@positions[mrk_index]
}

tidy_ehhs <- function(ehh_obj, pop_label) {
  # Works with recent rehh where as.data.frame(ehh_obj) returns POSITION, EHHS, NHAPLO, etc.
  df <- as.data.frame(ehh_obj)
  # If your column names differ (e.g., "position" or "ehhs"), rename below.
  stopifnot(all(c("POSITION") %in% names(df)))
  val_col <- if ("EHHS" %in% names(df)) "EHHS" else if ("ehhs" %in% names(df)) "ehhs" else stop("EHHS column not found")
  df %>%
    transmute(POSITION = .data$POSITION,
              value    = .data[[val_col]],
              POP      = pop_label,
              metric   = "EHHS")
}

tidy_ines <- function(scan_obj, pop_label) {
  df <- as.data.frame(scan_obj)
  # rehh's scan_hh() typically provides integrated site homozygosity as IES (aka inES)
  val_col <- intersect(c("IES","iES","IEHHS","INTEGRATED_EHHS"), names(df))
  stopifnot(length(val_col) == 1)
  pos_col <- intersect(c("POSITION","position","POS"), names(df))[1]
  df %>%
    transmute(POSITION = .data[[pos_col]],
              inES     = .data[[val_col]],
              POP      = pop_label)
}

tidy_rsb <- function(rsb_obj) {
  df <- as.data.frame(rsb_obj)
  # rehh::ies2rsb() returns columns like POSITION and Rsb (capitalization may vary)
  pos_col <- intersect(c("POSITION","position","POS"), names(df))[1]
  rsb_col <- intersect(c("Rsb","RSB","rsb"), names(df))[1]
  df %>%
    transmute(POSITION = .data[[pos_col]],
              Rsb      = .data[[rsb_col]])
}

# --- assemble `res` for your plotting function ---------------------------

make_res_for_plot <- function(hh_BP1, hh_TP1,
                              scan_BP1, scan_TP1,
                              ehh_BP1, ehh_TP1,
                              rsb_BP1_TP1,
                              mrk_index,
                              labels = c("BP1","TP1"),
                              chrom  = NULL) {
  # center coordinate from the SAME marker index in both subsets
  center <- get_center_pos(hh_BP1, mrk_index)
  
  # Panel A: EHHS curves as tidy long df
  ehhs_df <- bind_rows(
    tidy_ehhs(ehh_BP1, labels[1]),
    tidy_ehhs(ehh_TP1, labels[2])
  )
  
  # Panel B: per-pop inES
  ines_df <- bind_rows(
    tidy_ines(scan_BP1, labels[1]),
    tidy_ines(scan_TP1, labels[2])
  )
  
  # Panel C: Rsb (between pops)
  rsb_df <- tidy_rsb(rsb_BP1_TP1)
  
  # x-limits (window) — use the span of inES if you don't have a predefined window
  xmin <- min(ines_df$POSITION, na.rm = TRUE)
  xmax <- max(ines_df$POSITION, na.rm = TRUE)
  
  list(
    meta = list(
      chrom  = chrom %||% (if (!is.null(hh_BP1@chr.name)) hh_BP1@chr.name else NA_character_),
      idx_bp = center,
      window = c(xmin, xmax),
      labels = labels
    ),
    df = list(
      ehhs = ehhs_df,
      ines = ines_df,
      rsb  = rsb_df
    )
  )
}
