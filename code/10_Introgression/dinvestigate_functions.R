library(tidyverse)

# --- 95% empirical CI around the median
ci_quantile <- function(x, conf = 0.95) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(est = NA, lo = NA, hi = NA, n = 0))
  q <- quantile(x, probs = c((1 - conf)/2, 1 - (1 - conf)/2), na.rm = TRUE, names = FALSE)
  c(est = median(x, na.rm = TRUE), lo = q[1], hi = q[2], n = length(x))
}

# --- robust parser for "Comparision" (handles BP_ANC, TP_OUT)
parse_comp <- function(cmp) {
  toks <- strsplit(cmp, "_")[[1]]
  P1 <- toks[1]
  if (length(toks) >= 4 && paste(toks[(length(toks)-1):length(toks)], collapse = "_") == "BP_ANC") {
    P3 <- "BP_ANC"
    spacer <- paste(toks[2:(length(toks)-2)], collapse = "_")
  } else {
    P3 <- toks[length(toks)]
    spacer <- paste(toks[2:(length(toks)-1)], collapse = "_")
  }
  tibble(Comparision = cmp, P1 = P1, Spacer = spacer, P3 = P3)
}

# --- build "representative side" map: choose A-side comp (A = lower-index endpoint)
make_rep_map <- function(dinv_raw) {
  bp_order <- c("BP1","BP2","BP3","BP_ANC")
  comps <- tibble(Comparision = unique(dinv_raw$Comparision)) %>%
    rowwise() %>% do(parse_comp(.$Comparision)) %>% ungroup() %>%
    mutate(
      # canonical endpoints (A < B by bp_order)
      A = map2_chr(P1, P3, ~ c(.x, .y)[order(match(c(.x, .y), bp_order))][1]),
      B = map2_chr(P1, P3, ~ c(.x, .y)[order(match(c(.x, .y), bp_order))][2]),
      Hypothesis = paste0(A, "↔", B)
    )
  
  reps <- comps %>%
    group_by(Hypothesis) %>%
    # representative = row where P1 == A (the A-side comparison)
    arrange(desc(P1 == A)) %>% slice(1) %>%
    ungroup() %>%
    mutate(
      # trio label for inset
      trio_label = paste0("(", P1, ", ", Spacer, ")(", P3, ")"),
      # for A-side, f_dM<0 indicates A↔B sharing; flip sign so positive = A↔B
      sign_flip = -1
    ) %>%
    select(Hypothesis, Comparision, P1, Spacer, P3, trio_label, sign_flip)
  
  reps
}

# --- summarize using representative side only; flip sign so + = endpoint↔endpoint sharing
summarize_rep_side_globalbg <- function(dinv_raw, peaks.bed, rep_map, bg_stat = median) {
  peaks_tbl <- peaks.bed %>%
    transmute(chr = as.character(chrom), start, end,
              peak_num = as.numeric(peak), peaks = as.character(peak_num))
  
  out <- map_dfr(seq_len(nrow(rep_map)), function(i) {
    row <- rep_map[i, ]
    cmp <- row$Comparision
    flip <- row$sign_flip
    
    dcomp <- dinv_raw %>% filter(Comparision == cmp) %>%
      mutate(chr = as.character(chr),
             fd_adj = flip * f_dM)  # flip sign so + = A↔B sharing
    
    # genome-wide background on flipped values, excluding any peak windows
    if (nrow(dcomp) == 0) {
      bg_global <- NA_real_
    } else {
      ov_any <- rep(FALSE, nrow(dcomp))
      for (chr_i in unique(peaks_tbl$chr)) {
        idx <- which(dcomp$chr == chr_i); if (!length(idx)) next
        pk_chr <- peaks_tbl %>% filter(chr == chr_i)
        if (nrow(pk_chr)) {
          for (j in seq_len(nrow(pk_chr))) {
            ov_any[idx] <- ov_any[idx] |
              (dcomp$windowStart[idx] <= pk_chr$end[j] &
                 dcomp$windowEnd[idx]   >= pk_chr$start[j])
          }
        }
      }
      w_bg <- dcomp[!ov_any, , drop = FALSE]
      bg_global <- if (nrow(w_bg)) bg_stat(w_bg$fd_adj, na.rm = TRUE) else bg_stat(dcomp$fd_adj, na.rm = TRUE)
    }
    
    # per-peak median + 95% empirical CI (on flipped values)
    pk_rows <- pmap_dfr(peaks_tbl, function(chr, start, end, peak_num, peaks, ...) {
      w <- dcomp %>% filter(chr == !!chr)
      if (nrow(w) == 0) {
        tibble(chr, peak_num, peaks, df_hat = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, nwin = 0)
      } else {
        ov <- (w$windowStart <= end) & (w$windowEnd >= start)
        x  <- w$fd_adj[ov]; nwin <- sum(ov)
        if (nwin < 1) {
          tibble(chr, peak_num, peaks, df_hat = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, nwin = 0)
        } else {
          ci <- ci_quantile(x)
          tibble(chr, peak_num, peaks,
                 df_hat = as.numeric(ci["est"]),
                 ci_lo  = as.numeric(ci["lo"]),
                 ci_hi  = as.numeric(ci["hi"]),
                 nwin = nwin)
        }
      }
    })
    
    pk_rows %>%
      mutate(Hypothesis = row$Hypothesis,
             facet_label = row$Hypothesis,  # one facet per hypothesis
             bg_global = bg_global,
             trio_inset = row$trio_label)
  })
  
  out
}