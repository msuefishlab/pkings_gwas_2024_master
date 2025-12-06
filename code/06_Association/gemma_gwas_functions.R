library(ggrastr)

# functions inspired by erik enbody, with modifications

## USED##
plink.order <- function(plink.path, input.var) {
  df.plink <- fread(plink.path)

  # Parse chr & position from SNP like "chr1_993362_T_C"
  snp_parts <- str_split_fixed(df.plink$SNP, "_", 4)

  df.plink <- df.plink %>%
    mutate(
      CHR_from_snp = snp_parts[, 1],
      ps           = as.numeric(snp_parts[, 2]),
      numeric_chr  = as.numeric(str_replace_all(CHR_from_snp, "chr|ups", "")) # "chr1"→1
    ) %>%
    arrange(numeric_chr, ps) %>%
    mutate(
      CHR = factor(CHR_from_snp),
      log_p = -log10(.data[[input.var]]),
      chr_labels = CHR
    )

  # Row index shared across TESTs: 1,2,3,... in the order of unique (numeric_chr, ps)
  key <- paste(df.plink$numeric_chr, df.plink$ps, sep = ":")
  df.plink$row <- as.integer(factor(key, levels = unique(key)))

  # return plain data.frame to match your other helpers
  df.plink <- as.data.frame(df.plink)
  df.plink
}

gemma.order <- function(gemma.path, input.var) {
  df.gemma <- fread(gemma.path)
  df.gemma <- df.gemma %>% mutate(numeric_chr = as.numeric(str_replace_all(chr, "chr|ups", "")))
  df.gemma <- df.gemma %>% dplyr::arrange(numeric_chr, ps)
  df.gemma$chr <- factor(df.gemma$chr)
  df.gemma$row <- 1:nrow(df.gemma)
  df.gemma <- as.data.frame(df.gemma)
  df.gemma$log_p <- -log10(df.gemma[, input.var])
  df.gemma$chr_labels <- df.gemma$chr
  chr_breaks <- df.gemma %>%
    group_by(chr) %>%
    dplyr::summarise(chr_breaks = mean(row))
  df.gemma
}

## USED##

peak_list_permutation_plink <- function(input.df, input.var, thresholdA, num_snps, test) {
  autosome_peaks <- input.df %>%
    filter(TEST == test) %>%
    filter(!!sym(input.var) > thresholdA)

  # create granges object for individual SNPs that are above the threshold
  range.autosome_peaks <- GRanges(autosome_peaks$chr, IRanges(as.numeric(autosome_peaks$ps), as.numeric(autosome_peaks$ps)))

  # merge windows by some distance to other SNPs. defining it here by 75kb
  reduce.range.autosome_peaks <- GenomicRanges::reduce(range.autosome_peaks, min.gapwidth = 75000) # merge closest window

  # remove ranges that are made up of only less than num_snps
  # Initialize a vector to store the number of SNPs in each peak
  snps_in_peaks <- integer(length(reduce.range.autosome_peaks))

  # Count SNPs in each peak
  for (i in seq_along(reduce.range.autosome_peaks)) {
    overlaps <- findOverlaps(range.autosome_peaks, reduce.range.autosome_peaks[i])
    snps_in_peaks[i] <- length(unique(queryHits(overlaps)))
  }

  # Filter peaks based on the number of SNPs they contain
  valid_peaks <- snps_in_peaks >= num_snps
  reduce.range.autosome_peaks <- reduce.range.autosome_peaks[valid_peaks]

  # give each "peak" a name
  names(reduce.range.autosome_peaks) <- 1:length(reduce.range.autosome_peaks)

  datalist <- list()

  # create a dataframe where each SNP gets an peak ID
  for (i in (1:length(reduce.range.autosome_peaks))) {
    x <- autosome_peaks[queryHits(findOverlaps(range.autosome_peaks, reduce.range.autosome_peaks[i, ])), ]
    x$peak <- i
    datalist[[i]] <- x # add it to your list
  }
  peaks.named <- do.call(rbind, datalist)

  return(peaks.named)
}

peak_list_permutation <- function(input.df, input.var, thresholdA, num_snps) {
  autosome_peaks <- input.df %>% filter(!!sym(input.var) > thresholdA)

  # create granges object for individual SNPs that are above the threshold
  range.autosome_peaks <- GRanges(autosome_peaks$chr, IRanges(as.numeric(autosome_peaks$ps), as.numeric(autosome_peaks$ps)))

  # merge windows by some distance to other SNPs. defining it here by 75kb
  reduce.range.autosome_peaks <- GenomicRanges::reduce(range.autosome_peaks, min.gapwidth = 75000) # merge closest window

  # remove ranges that are made up of only less than num_snps
  # Initialize a vector to store the number of SNPs in each peak
  snps_in_peaks <- integer(length(reduce.range.autosome_peaks))

  # Count SNPs in each peak
  for (i in seq_along(reduce.range.autosome_peaks)) {
    overlaps <- findOverlaps(range.autosome_peaks, reduce.range.autosome_peaks[i])
    snps_in_peaks[i] <- length(unique(queryHits(overlaps)))
  }

  # Filter peaks based on the number of SNPs they contain
  valid_peaks <- snps_in_peaks >= num_snps
  reduce.range.autosome_peaks <- reduce.range.autosome_peaks[valid_peaks]

  # give each "peak" a name
  names(reduce.range.autosome_peaks) <- 1:length(reduce.range.autosome_peaks)

  datalist <- list()

  # create a dataframe where each SNP gets an peak ID
  for (i in (1:length(reduce.range.autosome_peaks))) {
    x <- autosome_peaks[queryHits(findOverlaps(range.autosome_peaks, reduce.range.autosome_peaks[i, ])), ]
    x$peak <- i
    datalist[[i]] <- x # add it to your list
  }
  peaks.named <- do.call(rbind, datalist)

  return(peaks.named)
}


plink.peak.bed <- function(peaks.named, log_p_threshold) {
  peaks.named %>%
    mutate(
      ps    = as.numeric(ps),
      log_p = as.numeric(log_p)
    ) %>%
    group_by(numeric_chr, peak) %>%
    summarise(
      start_bp = min(ps, na.rm = TRUE),
      end_bp = max(ps, na.rm = TRUE),
      region_len = max(ps, na.rm = TRUE) - min(ps, na.rm = TRUE), # avoid start/end refs
      num_snps = n(),
      mean_log_p = mean(log_p, na.rm = TRUE),
      max_log_p = max(log_p, na.rm = TRUE),
      num_snps_above_threshold = sum(log_p > log_p_threshold, na.rm = TRUE),
      most_significant_bp = ps[which.max(log_p)],
      .groups = "drop"
    )
}

### USED###
gemma.peak.bed <- function(peaks.named, log_p_threshold) {
  # Create a bed file that has a summary of each of these regions
  peak.bed <- peaks.named %>%
    group_by(numeric_chr, peak) %>%
    summarise(
      start = min(ps),
      end = max(ps),
      length = end - start,
      num.snps = n(),
      mean.log_p = mean(log_p),
      max.log_p = max(log_p),
      num.snps.above.threshold = sum(log_p > log_p_threshold),
      most_significant_bp = ps[which.max(log_p)], # Get the position of the most significant SNP
      .groups = "keep"
    )
  return(peak.bed)
}

### USED###
gemma.peak.ROW.bed <- function(peaks.named) {
  # create a bed file that has a summary of each of these regions
  peak.bed <- peaks.named %>%
    group_by(chr, peak) %>%
    summarise(
      start = min(row),
      end = max(row),
      length = end - start,
      num.snps = n(),
      mean.log_p = mean(log_p),
      max.log_p = max(log_p),
      .groups = "keep"
    )
  return(peak.bed)
}

# manc2.labels <- function(df.in, input.var){
#   # Assuming df.in is your input dataframe
#   chr_breaks <- df.in %>%
#     group_by(numeric_chr) %>%
#     dplyr::summarise(chr_breaks = mean(row))
#
#   chrom.colors <- data.frame(chr=unique(df.in$chr),
#                              color.num = rep(1:2,length(unique(df.in$chr)))) %>%
#     distinct(chr, .keep_all = T)
#
#   df.in2 <- df.in %>%
#     left_join(chrom.colors, by = "chr") %>%
#     mutate(color.num = as.factor(color.num))
#
#   df.in2 %>%
#     ggplot(aes_string(x = "row", y = input.var, col = "color.num",size=input.var)) + theme_bw() +
#     theme(legend.position="none",
#           #panel.border=element_blank(),
#           panel.border = element_blank(), axis.line.x = element_line(), axis.line.y = element_line(),
#           #axis.title.x=element_blank(),
#           #axis.text.x = element_text(angle = 45, color = "black"),
#           #axis.text.x = element_blank(),
#           panel.grid = element_blank(),
#           panel.background = element_blank(),
#           panel.grid.major.y=element_blank(),
#           panel.grid.minor.y=element_blank(),
#           axis.title.y = element_text(size=6),
#           axis.text = element_text(size=6),
#           axis.ticks.x=element_blank(),
#           axis.ticks.y=element_line(size=0.2)) +
#     ggrastr::geom_point_rast(shape = 20, stroke = 0.2, alpha = 0.25) +
#     scale_color_manual(values=rep(c("grey30","grey70"))) +
#     #scale_color_manual(values = c(rep_len(c("grey30", "red"), length(unique(chr_breaks$chr_ordered))+1))) #+
#     #scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1),minor_breaks = NULL) +
#     scale_x_continuous(
#       breaks = chr_breaks$chr_breaks,
#       labels = chr_breaks$numeric_chr
#     )+
#
#     labs(y=expression(paste("-log"[10], italic("P"),"-value")))
#   #scale_x_continuous(breaks=chr_breaks$chr_breaks,
#   #                   labels = chr_breaks$chr_labels)
# }

manc2.labels.plink <- function(df.in, input.var, chr.in = NULL, start = NULL, end = NULL) {
  library(dplyr)
  library(ggplot2)
  library(ggrastr)

  # 1) Subset by chromosome if requested
  if (!is.null(chr.in)) {
    df.in <- df.in %>% filter(CHR == chr.in)
  }

  # 2) Subset by window if both start and end are given
  if (!is.null(start) && !is.null(end)) {
    df.in <- df.in %>% filter(row >= start, row <= end)
  }

  # 3) Compute the x-axis breaks for whatever remains
  chr_breaks <- df.in %>%
    group_by(numeric_chr) %>%
    summarize(chr_breaks = mean(row), .groups = "drop")

  # 4) Recompute colors (will alternate grey30/grey70)
  chrom.colors <- tibble(
    chr = unique(df.in$chr),
    color.num = rep(1:2, length.out = length(unique(df.in$chr)))
  )

  df.in2 <- df.in %>%
    left_join(chrom.colors, by = c("CHR" = "chr")) %>%
    mutate(color.num = factor(color.num))

  # 5) Build the ggplot
  p <- ggplot(df.in2, aes_string(
    x    = "row",
    y    = input.var,
    col  = "color.num",
    size = input.var
  )) +
    theme_bw() +
    theme(
      legend.position   = "none",
      panel.border      = element_blank(),
      axis.line.x       = element_line(),
      axis.line.y       = element_line(),
      panel.grid        = element_blank(),
      axis.title.y      = element_text(size = 6),
      axis.text         = element_text(size = 6),
      axis.ticks.x      = element_blank(),
      axis.ticks.y      = element_line(size = 0.2)
    ) +
    ggrastr::geom_point_rast(
      shape = 20, stroke = 0.2, alpha = 0.25
    ) +
    scale_color_manual(
      values = rep(c("grey30", "grey70"), length.out = nrow(chr_breaks))
    ) +
    scale_x_continuous(
      breaks = chr_breaks$chr_breaks,
      labels = chr_breaks$numeric_chr
    ) +
    labs(
      y = expression(paste("-log"[10], italic("P"), "-value"))
    )

  return(p)
}

manc2.labels <- function(df.in, input.var, chr.in = NULL, start = NULL, end = NULL) {
  library(dplyr)
  library(ggplot2)
  library(ggrastr)

  # 1) Subset by chromosome if requested
  if (!is.null(chr.in)) {
    df.in <- df.in %>% filter(chr == chr.in)
  }

  # 2) Subset by window if both start and end are given
  if (!is.null(start) && !is.null(end)) {
    df.in <- df.in %>% filter(row >= start, row <= end)
  }

  # 3) Compute the x-axis breaks for whatever remains
  chr_breaks <- df.in %>%
    group_by(numeric_chr) %>%
    summarize(chr_breaks = mean(row), .groups = "drop")

  # 4) Recompute colors (will alternate grey30/grey70)
  chrom.colors <- tibble(
    chr = unique(df.in$chr),
    color.num = rep(1:2, length.out = length(unique(df.in$chr)))
  )

  df.in2 <- df.in %>%
    left_join(chrom.colors, by = "chr") %>%
    mutate(color.num = factor(color.num))

  # 5) Build the ggplot
  p <- ggplot(df.in2, aes_string(
    x    = "row",
    y    = input.var,
    col  = "color.num",
    size = input.var
  )) +
    theme_bw() +
    theme(
      legend.position   = "none",
      panel.border      = element_blank(),
      axis.line.x       = element_line(),
      axis.line.y       = element_line(),
      panel.grid        = element_blank(),
      axis.title.y      = element_text(size = 6),
      axis.text         = element_text(size = 6),
      axis.ticks.x      = element_blank(),
      axis.ticks.y      = element_line(size = 0.2)
    ) +
    ggrastr::geom_point_rast(
      shape = 20, stroke = 0.2, alpha = 0.25
    ) +
    scale_color_manual(
      values = rep(c("grey30", "grey70"), length.out = nrow(chr_breaks))
    ) +
    scale_x_continuous(
      breaks = chr_breaks$chr_breaks,
      labels = chr_breaks$numeric_chr
    ) +
    labs(
      y = expression(paste("-log"[10], italic("P"), "-value"))
    )

  return(p)
}

### CIRCOS-STYLE MANHATTAN PLOT ###
circos.manhattan <- function(df.in, input.var,
                             peaks_df = NULL,
                             sig_threshold = NULL,
                             sug_threshold = NULL,
                             plot_radius = 1.0,
                             gap_size = 1,
                             grid_lines = TRUE,
                             dynamic_max = NULL,
                             highlight_peaks = TRUE,
                             show_ideogram = TRUE,
                             # Color controls
                             chr_colors = c("grey30", "grey70"),
                             peak_color = "pink",
                             peak_alpha = 0.3,
                             peak_point_color = "red",
                             sig_line_color = "red",
                             sug_line_color = "blue",
                             grid_color = "grey90",
                             # Point size controls
                             point_size_range = c(0.1, 1.5),
                             point_alpha = 0.25,
                             peak_point_alpha = 0.4,
                             # Font size controls
                             chr_label_size = 0.8,
                             peak_label_size = 1.0,
                             axis_label_size = 0.7,
                             # Track height controls
                             ideogram_height = 0.05,
                             data_track_height = 0.7,
                             # Performance controls
                             downsample_low_p = TRUE,
                             downsample_threshold = 2,
                             downsample_frac = 0.1,
                             # Internal: for nested plot compatibility
                             .set_canvas_limits = TRUE) {
  library(circlize)
  library(dplyr)

  # Calculate dynamic max if not provided (99.999th percentile)
  if (is.null(dynamic_max)) {
    dynamic_max <- quantile(df.in[[input.var]], 0.99999, na.rm = TRUE)
  }

  y_limit <- ceiling(dynamic_max) + 2

  # Prepare chromosome information from full dataset
  chr_info <- df.in %>%
    group_by(chr, numeric_chr) %>%
    summarize(
      chr_len = max(ps, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(numeric_chr = as.numeric(numeric_chr)) %>% # Ensure numeric for proper sorting
    arrange(numeric_chr)

  # Prepare data for plotting
  plot_data <- df.in %>%
    mutate(log_p_value = .data[[input.var]]) %>%
    filter(!is.na(log_p_value), log_p_value >= 0, log_p_value <= y_limit)

  # Optional downsampling for performance
  if (downsample_low_p && nrow(plot_data) > 100000) {
    high_p <- plot_data %>% filter(log_p_value > downsample_threshold)
    low_p <- plot_data %>%
      filter(log_p_value <= downsample_threshold) %>%
      sample_frac(downsample_frac)
    plot_data <- bind_rows(high_p, low_p)
  }

  # Validate track heights don't exceed available space
  # Count tracks: data (1) + ideogram/labels (1)
  # When show_ideogram = FALSE, we create a 0.05 height label track
  n_tracks <- 2
  track_margin_total <- n_tracks * 0.01 # 0.005 on each side = 0.01 per track
  available_height <- 0.98 - track_margin_total

  # Label track is 0.05 when ideogram is hidden, otherwise use ideogram_height
  label_track_height <- if (show_ideogram) ideogram_height else 0.05
  total_track_height <- data_track_height + label_track_height

  if (total_track_height > available_height) {
    stop(paste0(
      "Total track height (", round(total_track_height, 2),
      ") exceeds available space (", round(available_height, 2), "). ",
      "Reduce data_track_height or ideogram_height."
    ))
  }

  # Initialize circos plot
  circos.clear()

  # plot_radius controls overall plot size/diameter
  # Larger values = larger plot, smaller values = smaller plot
  # Typical range: 0.5 (small) to 2 (large), default 1.0
  # Invert for canvas limits (larger plot_radius = smaller canvas limits)
  canvas_limit <- 1.0 / plot_radius

  # Set circos parameters
  # For nested plots, canvas limits should NOT be set
  if (.set_canvas_limits) {
    circos.par(
      start.degree = 90, # Start at top (12 o'clock)
      gap.degree = gap_size, # Gap between chromosomes
      track.margin = c(0.005, 0.005), # Smaller margin to fit tracks
      cell.padding = c(0, 0, 0, 0), # Padding within cells
      canvas.xlim = c(-canvas_limit, canvas_limit), # Plot diameter control
      canvas.ylim = c(-canvas_limit, canvas_limit),
      points.overflow.warning = FALSE # Suppress warnings for points outside region
    )
  } else {
    circos.par(
      start.degree = 90, # Start at top (12 o'clock)
      gap.degree = gap_size, # Gap between chromosomes
      track.margin = c(0.005, 0.005), # Smaller margin to fit tracks
      cell.padding = c(0, 0, 0, 0), # Padding within cells
      points.overflow.warning = FALSE # Suppress warnings for points outside region
    )
  }

  # Initialize with chromosome sectors
  # Set factor levels explicitly to maintain numeric order
  circos.initialize(
    factors = factor(chr_info$chr, levels = chr_info$chr),
    xlim = matrix(c(rep(0, nrow(chr_info)), chr_info$chr_len), ncol = 2)
  )

  # Track 1: GWAS points (outer track)
  circos.track(
    ylim = c(0, y_limit),
    track.height = data_track_height,
    bg.border = NA,
    panel.fun = function(x, y) {
      chr <- CELL_META$sector.index
      chr_num <- chr_info$numeric_chr[chr_info$chr == chr]

      # Get data for this chromosome
      chr_data <- plot_data %>% filter(chr == !!chr)

      if (nrow(chr_data) > 0) {
        # Point color based on chromosome
        color_idx <- ((chr_num - 1) %% length(chr_colors)) + 1
        point_col <- chr_colors[color_idx]

        # Point size based on log_p
        point_sizes <- scales::rescale(
          chr_data$log_p_value,
          to = point_size_range
        )

        # Plot points
        circos.points(
          chr_data$ps,
          chr_data$log_p_value,
          col = adjustcolor(point_col, alpha.f = point_alpha),
          pch = 16,
          cex = point_sizes
        )
      }
    }
  )

  # Track 2: Chromosome labels (when ideogram is hidden) or ideogram bands
  if (!show_ideogram) {
    # Create a thin track just for chromosome labels
    circos.track(
      ylim = c(0, 1),
      track.height = 0.05, # Small track just for labels
      bg.border = NA,
      panel.fun = function(x, y) {
        chr <- CELL_META$sector.index
        chr_num <- chr_info$numeric_chr[chr_info$chr == chr]

        # Add chromosome label centered in this track
        circos.text(
          CELL_META$xcenter,
          0.5,
          as.character(chr_num), # Convert to character for text display
          cex = chr_label_size,
          font = 2,
          facing = "bending.inside",
          niceFacing = TRUE
        )
      }
    )
  } else if (show_ideogram) {
    circos.track(
      ylim = c(0, 1),
      track.height = ideogram_height,
      bg.border = "white",
      panel.fun = function(x, y) {
        chr <- CELL_META$sector.index
        chr_num <- chr_info$numeric_chr[chr_info$chr == chr]

        # Alternating colors
        color_idx <- ((chr_num - 1) %% length(chr_colors)) + 1
        bg_color <- chr_colors[color_idx]

        # Fill the ideogram
        circos.rect(
          CELL_META$xlim[1], 0,
          CELL_META$xlim[2], 1,
          col = bg_color,
          border = "white",
          lwd = 0.5
        )

        # Add chromosome label (curved text)
        circos.text(
          CELL_META$xcenter,
          0.5,
          as.character(chr_num), # Convert to character for text display
          cex = chr_label_size,
          font = 2,
          facing = "bending.inside", # Text follows circle
          niceFacing = TRUE # Auto-flip for readability
        )
      }
    )
  }

  # Add grid lines and y-axis labels
  if (grid_lines) {
    grid_breaks <- seq(0, y_limit, by = 1)

    # Add y-axis to first chromosome only (cleaner look) - only on data track
    circos.yaxis(
      side = "left",
      sector.index = chr_info$chr[1],
      track.index = 1, # Only on the data track
      at = grid_breaks,
      labels = grid_breaks,
      labels.cex = axis_label_size,
      tick.length = 2,
      lwd = 0.3,
      labels.col = "black"
    )

    # Add circular grid lines
    for (grid_val in grid_breaks[-1]) { # Skip 0
      for (chr in chr_info$chr) {
        circos.lines(
          c(0, chr_info$chr_len[chr_info$chr == chr]),
          c(grid_val, grid_val),
          col = grid_color,
          lwd = 0.3,
          sector.index = chr,
          track.index = 1
        )
      }
    }
  }

  # Add significance threshold lines (circular rings)
  if (!is.null(sig_threshold)) {
    for (chr in chr_info$chr) {
      circos.lines(
        c(0, chr_info$chr_len[chr_info$chr == chr]),
        c(sig_threshold, sig_threshold),
        col = sig_line_color,
        lwd = 1.5,
        lty = 2,
        sector.index = chr,
        track.index = 1
      )
    }
  }

  if (!is.null(sug_threshold)) {
    for (chr in chr_info$chr) {
      circos.lines(
        c(0, chr_info$chr_len[chr_info$chr == chr]),
        c(sug_threshold, sug_threshold),
        col = sug_line_color,
        lwd = 1.5,
        lty = 2,
        sector.index = chr,
        track.index = 1
      )
    }
  }

  # Add peak highlighting
  if (!is.null(peaks_df) && highlight_peaks) {
    # Calculate peak boundaries
    peak_bounds <- peaks_df %>%
      group_by(peak, chr) %>%
      summarise(
        start_bp_raw = min(ps, na.rm = TRUE) - 50000, # Extend 50kb
        end_bp_raw = max(ps, na.rm = TRUE) + 50000,
        mid_bp = (min(ps, na.rm = TRUE) + max(ps, na.rm = TRUE)) / 2,
        .groups = "drop"
      ) %>%
      left_join(chr_info %>% select(chr, chr_len), by = "chr") %>%
      mutate(
        start_bp = pmax(start_bp_raw, 0),
        end_bp = pmin(end_bp_raw, chr_len)
      ) %>%
      select(-start_bp_raw, -end_bp_raw, -chr_len)

    # Add peak rectangles
    for (i in 1:nrow(peak_bounds)) {
      pb <- peak_bounds[i, ]

      circos.rect(
        pb$start_bp, 0,
        pb$end_bp, y_limit,
        col = adjustcolor(peak_color, alpha.f = peak_alpha),
        border = NA,
        sector.index = pb$chr,
        track.index = 1
      )

      # Add peak label
      circos.text(
        pb$mid_bp,
        y_limit * 0.95,
        pb$peak,
        cex = peak_label_size,
        font = 2,
        sector.index = pb$chr,
        track.index = 1,
        facing = "inside",
        niceFacing = TRUE
      )
    }

    # Highlight peak points in different color
    if (!is.null(sug_threshold)) {
      peak_points <- peaks_df %>%
        filter(
          .data[[input.var]] > (sug_threshold - 2),
          !is.na(.data[[input.var]]),
          .data[[input.var]] <= y_limit
        )

      for (chr in unique(peak_points$chr)) {
        chr_peaks <- peak_points %>% filter(chr == !!chr)

        if (nrow(chr_peaks) > 0) {
          point_sizes <- scales::rescale(
            chr_peaks[[input.var]],
            to = point_size_range
          )

          circos.points(
            chr_peaks$ps,
            chr_peaks[[input.var]],
            col = adjustcolor(peak_point_color, alpha.f = peak_point_alpha),
            pch = 16,
            cex = point_sizes,
            sector.index = chr,
            track.index = 1
          )
        }
      }
    }
  }

  # Return invisible NULL (plot is already drawn to device)
  invisible(NULL)
}

### CIRCOS-STYLE MANHATTAN PLOT WITH NESTED ZOOM ###
circos.manhattan.zoom <- function(df.in, input.var,
                                  zoom_peaks,
                                  zoom_extend = 0,
                                  parent_gap_degree = 1,
                                  zoom_gap_degree = 15,
                                  zoom_sector_width = 1,
                                  connection_col = "grey50",
                                  connection_height = 0.05,
                                  sig_threshold = NULL,
                                  sug_threshold = NULL,
                                  chr_colors = c("grey30", "grey70"),
                                  point_alpha = 0.25,
                                  parent_point_size = 0.3,
                                  zoom_point_size = 0.5,
                                  parent_label_size = 0.7,
                                  zoom_label_size = 0.6,
                                  zoom_xaxis_label_size = 0.5,
                                  parent_yaxis_label_size = 0.5,
                                  zoom_yaxis_label_size = 0.5,
                                  parent_ylim = NULL,
                                  zoom_ylim = NULL) {
  library(circlize)
  library(dplyr)

  # Capture input variable name
  input_var_str <- deparse(substitute(input.var))

  # Get chromosome info
  chr_info <- df.in %>%
    group_by(chr) %>%
    summarize(chr_len = max(ps), .groups = "drop") %>%
    mutate(numeric_chr = as.numeric(gsub("chr", "", chr))) %>%
    arrange(numeric_chr)

  # Create sector info for parent plot
  # Keep numeric chromosome order
  sector <- chr_info %>%
    arrange(numeric_chr) %>%
    select(chr, start = numeric_chr, end = chr_len) %>%
    mutate(start = 0) %>%
    select(chr, start, end) %>%
    as.data.frame()

  # Create correspondence table
  # Sort by chromosome number to match parent plot order
  correspondance <- zoom_peaks %>%
    mutate(parent_chr = paste0("chr", numeric_chr)) %>%
    left_join(chr_info %>% select(chr, chr_len),
              by = c("parent_chr" = "chr")) %>%
    mutate(
      parent_start = pmax(start - zoom_extend, 0),
      parent_end = pmin(end + zoom_extend, chr_len),
      zoom_name = paste0("Peak", peak, ":chr", numeric_chr),
      zoom_start = 0,
      zoom_end = zoom_sector_width  # Fixed width for all sectors
    ) %>%
    arrange(numeric_chr) %>%  # Sort by chromosome number
    select(parent_chr, parent_start, parent_end,
           zoom_name, zoom_start, zoom_end) %>%
    as.data.frame()

  # Create zoom sector info
  zoom_sector <- correspondance %>%
    select(zoom_name, zoom_start, zoom_end) %>%
    as.data.frame()

  # Filter data for zoom regions with coordinate scaling
  zoom_data_list <- list()
  for (i in 1:nrow(correspondance)) {
    # Get actual bp range for this region
    actual_start <- correspondance$parent_start[i]
    actual_end <- correspondance$parent_end[i]
    region_length <- actual_end - actual_start

    region_data <- df.in %>%
      filter(chr == correspondance$parent_chr[i],
             ps >= actual_start,
             ps <= actual_end) %>%
      mutate(
        name = correspondance$zoom_name[i],
        # Scale x-coordinates to fixed width
        x = (ps - actual_start) * zoom_sector_width / region_length,
        y = .data[[input_var_str]],
        # Store actual bp range for axis labeling
        actual_bp_start = actual_start,
        actual_bp_end = actual_end
      ) %>%
      select(name, x, y, actual_bp_start, actual_bp_end)
    zoom_data_list[[i]] <- region_data
  }
  zoom_data <- bind_rows(zoom_data_list) %>% as.data.frame()

  # Parent data
  parent_data <- df.in %>%
    select(chr, x = ps, y = all_of(input_var_str)) %>%
    as.data.frame()

  # Determine y-axis limits with padding
  # If not provided, calculate automatically
  if (is.null(parent_ylim)) {
    parent_ylim <- c(0, ceiling(max(df.in[[input_var_str]], na.rm = TRUE)) + 2)
  }

  if (is.null(zoom_ylim)) {
    zoom_ylim <- c(0, ceiling(max(df.in[[input_var_str]], na.rm = TRUE)) + 2)
  }

  # Parent plot function
  f1 <- function() {
    circos.par(gap.degree = parent_gap_degree, start.degree = 90, points.overflow.warning = FALSE)
    circos.initialize(factors = factor(sector$chr, levels = sector$chr),
                     xlim = sector[, 2:3])
    circos.track(parent_data[[1]],
                x = parent_data[[2]],
                y = parent_data[[3]],
                ylim = parent_ylim,
                bg.border = NA,
                panel.fun = function(x, y) {
                  # Get chromosome number for color alternation
                  chr <- CELL_META$sector.index
                  chr_num <- chr_info$numeric_chr[chr_info$chr == chr]
                  color_idx <- ((chr_num - 1) %% length(chr_colors)) + 1
                  base_color <- chr_colors[color_idx]

                  # Plot points
                  point_colors <- ifelse(
                    !is.null(sig_threshold) & y > sig_threshold,
                    "red",
                    base_color
                  )
                  circos.points(x, y,
                               col = adjustcolor(point_colors,
                                               alpha.f = point_alpha),
                               pch = 16, cex = parent_point_size)

                  # Threshold lines
                  if (!is.null(sig_threshold)) {
                    circos.lines(CELL_META$xlim,
                               c(sig_threshold, sig_threshold),
                               col = "red", lwd = 2, lty = 2)
                  }
                  if (!is.null(sug_threshold)) {
                    circos.lines(CELL_META$xlim,
                               c(sug_threshold, sug_threshold),
                               col = "blue", lwd = 2, lty = 2)
                  }

                  # Labels - outside the circle at top of track
                  circos.text(CELL_META$xcenter,
                            CELL_META$ylim[2] + mm_y(2),
                            CELL_META$sector.index,
                            niceFacing = TRUE,
                            adj = c(0.5, 0), cex = parent_label_size)
                })

    # Add y-axis to first chromosome only
    circos.yaxis(side = "left",
                sector.index = sector$chr[1],
                labels.cex = parent_yaxis_label_size)

    # Add horizontal grid lines across all chromosomes
    grid_breaks <- seq(parent_ylim[1], parent_ylim[2], by = 1)
    for (grid_val in grid_breaks[-1]) {  # Skip the first value (y=0)
      for (chr in sector$chr) {
        circos.lines(
          c(0, chr_info$chr_len[chr_info$chr == chr]),
          c(grid_val, grid_val),
          col = "grey90",
          lwd = 0.3,
          sector.index = chr,
          track.index = 1
        )
      }
    }
  }

  # Zoom plot function
  f2 <- function() {
    circos.par(gap.degree = zoom_gap_degree,
              cell.padding = c(0, 0, 0, 0),
              points.overflow.warning = FALSE)
    circos.initialize(
      factors = factor(zoom_sector$zoom_name, levels = zoom_sector$zoom_name),
      xlim = as.matrix(zoom_sector[, 2:3]))
    circos.track(zoom_data[[1]],
                x = zoom_data[[2]],
                y = zoom_data[[3]],
                ylim = zoom_ylim,
                panel.fun = function(x, y) {
                  point_colors <- ifelse(
                    !is.null(sig_threshold) & y > sig_threshold,
                    "red",
                    chr_colors[1]
                  )
                  circos.points(x, y,
                               col = adjustcolor(point_colors,
                                               alpha.f = point_alpha),
                               pch = 16, cex = zoom_point_size)

                  # Threshold lines
                  if (!is.null(sig_threshold)) {
                    circos.lines(CELL_META$xlim,
                               c(sig_threshold, sig_threshold),
                               col = "red", lwd = 0.5, lty = 2)
                  }
                  if (!is.null(sug_threshold)) {
                    circos.lines(CELL_META$xlim,
                               c(sug_threshold, sug_threshold),
                               col = "blue", lwd = 0.5, lty = 2)
                  }

                  # Add x-axis with Mb-scaled positions - beneath the plot
                  # Get actual bp range for this sector
                  if (length(x) > 0) {
                    sector_zoom_data <- zoom_data[zoom_data$name ==
                                                  CELL_META$sector.index, ]
                    if (nrow(sector_zoom_data) > 0) {
                      actual_start <- sector_zoom_data$actual_bp_start[1]
                      actual_end <- sector_zoom_data$actual_bp_end[1]

                      # Generate nice Mb tick positions
                      mb_ticks <- pretty(c(actual_start, actual_end) / 1e6,
                                        n = 5)
                      mb_ticks_bp <- mb_ticks * 1e6

                      # Scale to fixed coordinate system
                      region_length <- actual_end - actual_start
                      scaled_ticks <- (mb_ticks_bp - actual_start) *
                                     zoom_sector_width / region_length

                      # Filter to ticks within range
                      valid_idx <- scaled_ticks >= 0 &
                                  scaled_ticks <= zoom_sector_width
                      scaled_ticks <- scaled_ticks[valid_idx]
                      mb_labels <- paste0(mb_ticks[valid_idx], "M")

                      # Create axis with Mb labels - on outer edge
                      circos.axis(h = "bottom",
                                 major.at = scaled_ticks,
                                 labels = mb_labels,
                                 labels.cex = zoom_xaxis_label_size,
                                 labels.facing = "reverse.clockwise",
                                 direction="inside",
                                 minor.ticks = 0)
                    }
                  }

                  # Add peak label in interior, below x-axis
                  # Extract peak number from sector name (e.g., "Peak1:chr8" -> "Peak 1")
                  peak_num <- gsub("Peak([0-9]+):.*", "Peak \\1", CELL_META$sector.index)
                  circos.text(CELL_META$xcenter,
                            CELL_META$ylim[1] - mm_y(2),
                            peak_num,
                            facing = "bending.inside",
                            adj = c(0.5, 5), cex = zoom_label_size)
                },
                track.margin = c(0, 0))

    # Add y-axis to first zoom sector only
    circos.yaxis(side = "left",
                sector.index = zoom_sector$zoom_name[1],
                labels.cex = zoom_yaxis_label_size)

    # Add horizontal grid lines across all zoom sectors
    grid_breaks <- seq(zoom_ylim[1], zoom_ylim[2], by = 1)
    for (grid_val in grid_breaks[-1]) {  # Skip the first value (y=0)
      for (zoom_name in zoom_sector$zoom_name) {
        circos.lines(
          c(0, zoom_sector_width),
          c(grid_val, grid_val),
          col = "grey90",
          lwd = 0.3,
          sector.index = zoom_name,
          track.index = 1
        )
      }
    }
  }

  # Create nested plot
  # f1 (parent) outside, f2 (zoom) inside
  # connection_height: height of connection track (lower = closer circles)
  circos.nested(f1, f2, correspondance,
                connection_col = connection_col,
                connection_height = connection_height)
}
