
get_thr <- function(test, kind = c("suggestive","significant")){
  kind <- match.arg(kind)
  thr_by_test_ord %>% dplyr::filter(TEST == test) %>% dplyr::pull({{kind}})
}

safe_write_table <- function(df, path, col_names = TRUE){
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = col_names)
}

# core peak-buiding for a TEST
build_peaks_for_test <- function(df_all, test){
  # Filter genome-wide table to reduce memory footprint
  df <- df_all %>%
    dplyr::filter(log_p > 1,
                  MAF > params$MAF_thresh,
                  !is.na(log_p),
                  !grepl("ups", CHR)) %>%
    dplyr::arrange(numeric_chr)
  
  # SNP lists & peak BED (unfiltered)
  peaks_tbl <- peak_list_permutation_plink(df %>% dplyr::filter(TEST == test), 
                                           "log_p", 3, 4, test)  # (win params as before)
  bed_unfiltered <- plink.peak.bed(peaks_tbl, get_thr(test, "suggestive"))
  
  list(df=df, peaks_tbl=peaks_tbl, bed_unfiltered=bed_unfiltered)
}

# filter/renumber peaks & derive per-peak SNPs/topN/index
finalize_peaks <- function(peaks_tbl, bed_unfiltered, test, topN=250){
  thr_sugg <- get_thr(test, "suggestive")
  
  bed_f <- bed_unfiltered %>%
    dplyr::filter(max_log_p > thr_sugg,
                  num_snps_above_threshold > 0) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(numeric_chr, start_bp) %>%
    dplyr::mutate(
      peak.original = peak,
      peak.new = dplyr::row_number(),
      peak = peak.new
    )
  
  peak_name_map <- bed_f %>% dplyr::select(peak.new, peak.original)
  
  peaks_in <- peaks_tbl %>%
    dplyr::filter(peak %in% bed_f$peak.original) %>%
    dplyr::left_join(peak_name_map, by = c("peak" = "peak.original")) %>%
    dplyr::rename(peak.original = peak, peak = peak.new) %>%
    dplyr::arrange(peak)
  
  peaks_topN <- peaks_in %>%
    dplyr::group_by(peak) %>%
    dplyr::slice_max(log_p, n = topN, with_ties = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::select(CHR, ps, log_p, peak)
  
  index_snps <- peaks_in %>%
    dplyr::group_by(peak) %>%
    dplyr::filter(log_p == max(log_p)) %>%
    dplyr::slice_sample(n = 1) %>%
    dplyr::ungroup()
  
  list(bed_filtered = bed_f,
       peaks_in = peaks_in,
       peaks_topN = peaks_topN,
       index_snps = index_snps,
       thr_sugg = thr_sugg)
}

# write all artifacts for a TEST
write_artifacts <- function(test, out, outdir){
  tag <- tolower(test)
  # SNPs in peaks
  safe_write_table(out$peaks_in,
                   file.path(outdir, sprintf("%s_PLINK_%s_SNPS_IN_PEAKS.txt", params$data_path, test)))
  # TopN per peak
  safe_write_table(out$peaks_topN,
                   file.path(outdir, sprintf("%s_PLINK_%s_SNPS_IN_PEAKS_TOP250.txt", params$data_path, test)),
                   col_names = FALSE)
  # Index SNPs (one per peak)
  safe_write_table(out$index_snps$SNP,
                   file.path(outdir, sprintf("%s_PLINK_%s_index_snps.txt", params$data_path, test)),
                   col_names = FALSE)
  # BED unfiltered
  safe_write_table(out$bed_unfiltered %>% dplyr::select(numeric_chr, start_bp, end_bp, peak),
                   file.path(outdir, sprintf("%s_PLINK_%s_PEAKS_unfiltered.bed", params$data_path, test)),
                   col_names = FALSE)
  # BED filtered (ordered by max log_p for convenience)
  safe_write_table(out$bed_filtered %>% dplyr::arrange(max_log_p) %>%
                     dplyr::select(numeric_chr,start_bp,end_bp,peak,max_log_p,most_significant_bp),
                   file.path(outdir, sprintf("%s_PLINK_%s_PEAKS.bed", params$data_path, test)),
                   col_names = FALSE)
}

# Manhattan plot with shaded peaks (per TEST)
plot_manhattan_with_peaks <- function(df_all, peaks_in, test, thr_sugg, thr_sig){
  df_t <- df_all %>% dplyr::filter(TEST == test)
  dynamic_max <- stats::quantile(df_t$log_p, 0.99999, na.rm = TRUE)
  
  peak_bounds <- peaks_in %>%
    dplyr::group_by(peak) %>%
    dplyr::summarise(
      xmin = min(row, na.rm=TRUE) - 1e4,
      xmax = max(row, na.rm=TRUE) + 1e4,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      ymin = 1,
      ymax = dynamic_max
    )
  
  p <- manc2.labels.plink(df_t, "log_p") +
    ggplot2::geom_rect(
      data = peak_bounds,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "grey", alpha = 0.3, inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = peak_bounds,
      ggplot2::aes(x = (xmin + xmax)/2, y = dynamic_max + 1, label = peak),
      inherit.aes = FALSE, size = 4, color = "black"
    ) +
    ggplot2::geom_hline(yintercept = thr_sugg, linetype = 2, color = "blue") +
    ggplot2::geom_hline(yintercept = thr_sig,  linetype = 2, color = "red") +
    ggplot2::geom_point(
      data = dplyr::filter(peaks_in, log_p > thr_sugg),
      ggplot2::aes(x = row, y = log_p, size = log_p),
      color = "red", shape = 20, stroke = 0.2, alpha = 0.25, inherit.aes = FALSE
    ) +
    ggplot2::scale_size_continuous(range = c(0.9, 3)) +
    ggplot2::labs(y = expression(paste("-log"[10], italic("P"), "-value")), x = "Chromosome") +
    ggplot2::ylim(c(1, dynamic_max + 2)) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 14),
      axis.title.x = ggplot2::element_text(size = 16),
      axis.title.y = ggplot2::element_text(size = 16)
    )
  
  p
}

# End-to-end for one TEST
process_test <- function(test){
  built <- build_peaks_for_test(pk.lmm, test)
  fin   <- finalize_peaks(built$peaks_tbl, built$bed_unfiltered, test, topN = 250)
  fin$bed_unfiltered <- built$bed_unfiltered  # carry for writing
  write_artifacts(test, fin, peak_data_path)
  
  # Optional quick histogram of #SNPs/peak for QC
  print(
    ggplot2::ggplot(fin$bed_filtered, ggplot2::aes(x = num_snps)) +
      ggplot2::geom_histogram(bins = 10) + ggplot2::scale_x_log10() +
      ggplot2::ggtitle(sprintf("%s: SNPs per filtered peak", test))
  )
  
  # Manhattan
  p <- plot_manhattan_with_peaks(
    df_all   = pk.lmm,
    peaks_in = fin$peaks_in,
    test     = test,
    thr_sugg = get_thr(test, "suggestive"),
    thr_sig  = get_thr(test, "significant")
  )
  print(p)
  
  # Save SVG
  out_svg <- file.path(root,"output_data","06_Association", params$data_path,
                       sprintf("%s.gwas.plink.%s.svg", params$data_path, tolower(test)))
  ggsave(out_svg, plot = p, width = 12, height = 4, units = "in")
  
  invisible(list(test=test, results=fin, plot=p, svg=out_svg))
}