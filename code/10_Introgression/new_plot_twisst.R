library(ggplot2)
library(dplyr)
library(stringr)
library(ape)
library(purrr)
library(cowplot)
library(grid) # for grid::grobTree


## import_twisst_all: Refactored with per-entry embedded-topology loading and correct closure
# Reads topologies from an external Newick file if provided; otherwise extracts per-entry
import_twisst_all <- function(files_df,
                              topos_file = NULL,
                              min_subtrees = 1,
                              max_interval = Inf,
                              recalc_mid = TRUE) {
  # Validate input
  required_cols <- c("chrom", "comparison", "topocounts", "intervals")
  if (!all(required_cols %in% names(files_df))) {
    stop("files_df must have columns: ", paste(required_cols, collapse = ", "))
  }

  # Helper to read intervals & topocounts
  read_pair <- function(tc_file, iv_file) {
    iv <- read.table(iv_file, header = TRUE)
    tc <- read.table(tc_file, header = TRUE)
    list(
      interval_data = iv,
      topocounts    = tc,
      tc_file       = tc_file,
      iv_file       = iv_file
    )
  }

  # Import raw data
  df_imp <- files_df %>%
    mutate(
      topocounts = as.character(topocounts),
      intervals  = as.character(intervals)
    ) %>%
    rowwise() %>%
    mutate(raw = list(read_pair(topocounts, intervals))) %>%
    ungroup()

  # Process each dataset entry
  out_list <- df_imp %>%
    mutate(key = paste0(chrom, "__", comparison)) %>%
    select(key, raw) %>%
    deframe() %>%
    map(function(entry) {
      iv <- entry$interval_data
      tc <- entry$topocounts
      # Filter windows
      keep <- which(
        !is.na(rowSums(tc)) &
          rowSums(tc) >= min_subtrees &
          (iv$end - iv$start + 1) <= max_interval
      )
      iv <- iv[keep, , drop = FALSE]
      tc <- tc[keep, , drop = FALSE]

      # Compute weights and midpoints
      w <- tc / rowSums(tc)
      if (!"mid" %in% names(iv) || recalc_mid) iv$mid <- (iv$start + iv$end) / 2

      # Load per-entry embedded topologies or external
      topos_entry <- NULL
      if (is.null(topos_file)) {
        raw_tc <- read.table(entry$tc_file, header = TRUE)
        n_topos <- ncol(raw_tc) - 1
        newick_txt <- read.table(entry$tc_file,
          nrow = n_topos,
          comment.char = "",
          sep = "	",
          as.is = TRUE
        )[, 1]
        topos_entry <- tryCatch(read.tree(text = newick_txt), error = function(e) NULL)
      } else {
        topos_entry <- tryCatch(read.tree(file = topos_file), error = function(e) NULL)
      }

      # Build result list
      res <- list(
        interval_data = iv,
        topocounts = tc,
        weights = w,
        weights_mean = colMeans(w, na.rm = TRUE),
        weights_overall = colSums(w) / sum(colSums(w)),
        pos = iv$mid
      )
      if (!is.null(topos_entry)) res$topologies <- topos_entry
      res
    })

  return(out_list)
}

# Example usage:
# files_tbl <- tibble::tibble(...)
# twisst_data <- import_twisst_all(files_tbl, topos_file = NULL)  # or path to tree file
# files_tbl <- tibble::tibble(...)
# twisst_data <- import_twisst_all(files_tbl, topos_file=NULL)  # or path to tree file

# Helper to find indices of topologies where grp1+grp2 are monophyletic
find_keep_idx <- function(topos, grp1, grp2) {
  which(sapply(topos, function(tr) {
    labs <- tr$tip.label
    # ensure both groups exist in this tree
    if (!(grp1 %in% labs && grp2 %in% labs)) {
      return(FALSE)
    }
    ape::is.monophyletic(tr, c(grp1, grp2))
  }))
}

plot_twisst_region <- function(twisst_data, chrom, start, end,
                               tree_width = 1,
                               tree_height = 1) {
  # Define comparison colors
  colors <- c(
    bp1_bp2 = "#1b9e77",
    bp1_bp3 = "#d95f02",
    bp2_bp3 = "#7570b3"
  )

  # 1. Extract relevant keys for this chromosome and comparison types
  keys <- grep(paste0("^", chrom, "__.*(bp1_bp2|bp1_bp3|bp2_bp3)"),
    names(twisst_data),
    value = TRUE
  )

  # 2. Summarize data & extract one tree per comparison
  summary_list <- lapply(keys, function(k) {
    cmp <- sub(paste0(chrom, "__"), "", k)
    parts <- str_split(cmp, "_")[[1]]
    g1 <- toupper(parts[1])
    g2 <- toupper(parts[2])
    tb <- twisst_data[[k]]
    keep <- which(sapply(tb$topologies, function(tr) {
      labs <- tr$tip.label
      (g1 %in% labs && g2 %in% labs) &&
        ape::is.monophyletic(tr, c(g1, g2))
    }))
    if (length(keep) == 0) stop("No matching topology for ", cmp)
    wvec <- rowSums(tb$weights[, keep, drop = FALSE])
    list(
      comparison = cmp,
      data       = data.frame(pos = tb$pos, weight = wvec),
      tree       = tb$topologies[[keep[1]]]
    )
  })

  # 3. Combine and filter weight data for plotting
  df_all <- bind_rows(lapply(summary_list, function(x) {
    cbind(comparison = x$comparison, x$data)
  })) %>%
    filter(pos >= start, pos <= end)

  # 4. Helper: compute panel width based on number of edges
  get_panel_width <- function(tr) {
    (max(node.depth(tr)) - 1) * tree_width
  }

  # 5. Set up dynamic layout: top titles | middle trees | bottom weight plot
  n <- length(summary_list)
  layout(
    rbind(
      seq_len(n),
      seq(n + 1, 2 * n),
      rep(2 * n + 1, n)
    ),
    heights = c(0.2, tree_height, 1),
    widths = rep(1, n)
  )

  # Row 1: Titles
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  for (i in seq_len(n)) {
    plot.new()
    cmp <- summary_list[[i]]$comparison
    txt <- toupper(gsub("_", " vs ", cmp))
    text(0.5, 0.5, txt, col = colors[cmp], cex = 1)
  }

  # Row 2: Trees
  for (i in seq_len(n)) {
    tr <- summary_list[[i]]$tree
    cmp <- summary_list[[i]]$comparison
    # Vertical limits based on node heights
    nh <- node.height(tr) * tree_height
    ylim_tree <- c(min(nh) - 1.5 * tree_height * 0.1, max(nh) + 0.1 * tree_height)
    # Horizontal limits based on edge count
    panel_w <- get_panel_width(tr)

    plot(NA,
      xlim = c(0, panel_w),
      ylim = ylim_tree,
      xaxt = "n", yaxt = "n",
      xlab = "", ylab = "", bty = "n"
    )

    # Center and draw tree
    xmid <- panel_w / 2
    draw.tree(
      phy = tr,
      x = xmid,
      y = 0,
      x_scale = tree_width,
      y_scale = tree_height,
      col = colors[cmp],
      col.label = colors[cmp],
      add_labels = TRUE,
      label_offset = 0.1 * tree_height,
      cex = 1,
      lwd = 2,
      direction = "right"
    )
  }

  # Row 3: Weight plot
  par(mar = c(4, 4, 1, 1))
  first_cmp <- summary_list[[1]]$comparison
  plot(df_all$pos[df_all$comparison == first_cmp],
    df_all$weight[df_all$comparison == first_cmp],
    type = "l", col = colors[first_cmp],
    xlab = "Genomic position",
    ylab = "Weight",
    ylim = c(0, max(df_all$weight)),
    main = paste0(chrom, ": ", start, "-", end)
  )
  for (i in 2:n) {
    cmp <- summary_list[[i]]$comparison
    lines(df_all$pos[df_all$comparison == cmp],
      df_all$weight[df_all$comparison == cmp],
      col = colors[cmp], lwd = 1
    )
  }
  abline(h = 1 / 3, lty = "dashed", col = "gray50")

  # Reset layout
  layout(1)
}

plot_twisst_region_with_manhattan <- function(twisst_data, chrom, start, end,
                                              tree_width = 1,
                                              tree_height = 1,
                                              manhattan_df,
                                              manhattan_var,
                                              dynamic_max) {
  # Define comparison colors
  colors <- c(
    bp1_bp2 = "#1b9e77",
    bp1_bp3 = "#d95f02",
    bp2_bp3 = "#7570b3"
  )

  # 1. Extract relevant keys for this chromosome and comparison types
  keys <- grep(paste0("^", chrom, "__.*(bp1_bp2|bp1_bp3|bp2_bp3)"),
    names(twisst_data),
    value = TRUE
  )

  # 2. Summarize data & extract one tree per comparison
  summary_list <- lapply(keys, function(k) {
    cmp <- sub(paste0(chrom, "__"), "", k)
    parts <- str_split(cmp, "_")[[1]]
    g1 <- toupper(parts[1])
    g2 <- toupper(parts[2])
    tb <- twisst_data[[k]]
    keep <- which(sapply(tb$topologies, function(tr) {
      labs <- tr$tip.label
      (g1 %in% labs && g2 %in% labs) &&
        ape::is.monophyletic(tr, c(g1, g2))
    }))
    if (length(keep) == 0) stop("No matching topology for ", cmp)
    wvec <- rowSums(tb$weights[, keep, drop = FALSE])
    list(
      comparison = cmp,
      data       = data.frame(pos = tb$pos, weight = wvec),
      tree       = tb$topologies[[keep[1]]]
    )
  })

  # 3. Combine and filter weight data for plotting
  df_all <- bind_rows(lapply(summary_list, function(x) {
    cbind(comparison = x$comparison, x$data)
  })) %>%
    filter(pos >= start, pos <= end)

  # 4. Helper: compute panel width based on number of edges
  get_panel_width <- function(tr) {
    (max(node.depth(tr)) - 1) * tree_width
  }

  # 5. Set up dynamic layout: top titles | middle trees | weight plot | Manhattan plot
  n <- length(summary_list)
  layout(
    rbind(
      seq_len(n), # Titles
      seq(n + 1, 2 * n), # Trees
      rep(2 * n + 1, n), # Weight plot
      rep(2 * n + 2, n) # Manhattan plot
    ),
    heights = c(0.2, tree_height, 1, 1),
    widths = rep(1, n)
  )

  # Row 1: Titles
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  for (i in seq_len(n)) {
    plot.new()
    cmp <- summary_list[[i]]$comparison
    txt <- toupper(gsub("_", " vs ", cmp))
    text(0.5, 0.5, txt, col = colors[cmp], cex = 1)
  }

  # Row 2: Trees
  for (i in seq_len(n)) {
    tr <- summary_list[[i]]$tree
    cmp <- summary_list[[i]]$comparison
    # Vertical limits based on node heights
    nh <- node.height(tr) * tree_height
    ylim_tree <- c(min(nh) - 1.5 * tree_height * 0.1, max(nh) + 0.1 * tree_height)
    # Horizontal limits based on edge count
    panel_w <- get_panel_width(tr)

    plot(NA,
      xlim = c(0, panel_w),
      ylim = ylim_tree,
      xaxt = "n", yaxt = "n",
      xlab = "", ylab = "", bty = "n"
    )

    # Center and draw tree
    xmid <- panel_w / 2
    draw.tree(
      phy = tr,
      x = xmid,
      y = 0,
      x_scale = tree_width,
      y_scale = tree_height,
      col = colors[cmp],
      col.label = colors[cmp],
      add_labels = TRUE,
      label_offset = 0.1 * tree_height,
      cex = 1,
      lwd = 2,
      direction = "right"
    )
  }

  # Row 3: Weight plot
  par(mar = c(4, 4, 1, 1))
  first_cmp <- summary_list[[1]]$comparison
  plot(df_all$pos[df_all$comparison == first_cmp],
    df_all$weight[df_all$comparison == first_cmp],
    type = "l", col = colors[first_cmp],
    xlab = "Genomic position",
    ylab = "Weight",
    ylim = c(0, max(df_all$weight)),
    main = paste0(chrom, ": ", start, "-", end)
  )
  for (i in 2:n) {
    cmp <- summary_list[[i]]$comparison
    lines(df_all$pos[df_all$comparison == cmp],
      df_all$weight[df_all$comparison == cmp],
      col = colors[cmp], lwd = 1
    )
  }
  abline(h = 1 / 3, lty = "dashed", col = "gray50")

  par(mar = c(4, 4, 1, 1))
  df_sub <- manhattan_df
  if (!is.null(chrom)) {
    df_sub <- df_sub %>% filter(chr == chrom)
  }
  df_sub <- df_sub %>% filter(row >= start, row <= end)

  with(df_sub, {
    plot(row, get(manhattan_var),
      col = rgb(0.2, 0.2, 0.2, alpha = 0.25),
      pch = 16,
      cex = 1,
      xlab = "Chromosome",
      ylab = expression(paste("-log"[10], italic("P"), "-value")),
      main = "Manhattan Plot",
      ylim = c(1, dynamic_max + 2),
      axes = TRUE
    )
  })

  # Add grey rectangle for peak_bounds[2, ]
  rect(
    xleft = peak_bounds[2, "xmin"],
    xright = peak_bounds[2, "xmax"],
    ybottom = peak_bounds[2, "ymin"],
    ytop = dynamic_max,
    col = rgb(0.5, 0.5, 0.5, alpha = 0.3),
    border = NA
  )

  # Add label at top of peak
  text(
    x = (peak_bounds[2, "xmin"] + peak_bounds[2, "xmax"]) / 2,
    y = dynamic_max + 1,
    labels = peak_bounds[2, "peak"],
    cex = 1.2
  )

  # Add horizontal lines for thresholds
  abline(h = suggestive, col = "blue", lty = 2)
  abline(h = significant, col = "red", lty = 2)

  # Reset layout
  layout(1)
}

# Usage example:
# plot_twisst_region(my_twisst_data, "Chr1", 100000, 200000,
#                    tree_width = 2, tree_height = 1)


# function for plotting tree that uses ape to get node positions
draw.tree <- function(phy, x, y, x_scale = 1, y_scale = 1, method = 1, direction = "right",
                      col = "black", col.label = "black", add_labels = TRUE, add_symbols = FALSE,
                      label_offset = .02, symbol_offset = 0, col.symbol = "black", symbol_bg = "NA",
                      pch = 19, cex = NULL, lwd = NULL, label_alias = NULL) {
  n_tips <- length(phy$tip.label)

  if (direction == "right") {
    node_x <- (max(node.depth(phy, method = method)) - node.depth(phy, method = method))
    node_y <- node.height(phy) * y_scale
    label_x <- node_x[1:n_tips] + label_offset
    label_y <- node_y[1:n_tips]
    adj_x <- 0
    adj_y <- .5
    symbol_x <- node_x[1:n_tips] + symbol_offset
    symbol_y <- node_y[1:n_tips]
  }
  if (direction == "down") {
    node_y <- (node.depth(phy, method = method) - 1) * y_scale * 1
    node_x <- node.height(phy) * x_scale
    label_x <- node_x[1:n_tips]
    label_y <- node_y[1:n_tips] - label_offset
    adj_x <- .5
    adj_y <- 1
    symbol_x <- node_x[1:n_tips]
    symbol_y <- node_y[1:n_tips] - symbol_offset
  }

  # draw edges
  segments(x + node_x[phy$edge[, 1]], y + node_y[phy$edge[, 1]],
    x + node_x[phy$edge[, 2]], y + node_y[phy$edge[, 2]],
    col = col, lwd = lwd
  )

  if (is.null(label_alias) == FALSE) {
    tip_labels <- label_alias[phy$tip.label]
  } else {
    tip_labels <- phy$tip.label
  }

  if (add_labels == "TRUE") text(x + label_x, y + label_y, col = col.label, labels = tip_labels, adj = c(adj_x, adj_y), cex = cex)
  if (add_symbols == "TRUE") poin
  ts(x + symbol_x, y + symbol_y, pch = pch, col = col.symbol, bg = symbol_bg)
}


plot_twisst_region_gg <- function(twisst_data, chrom, start, end,
                                  tree_width = 1,
                                  tree_height = 1) {
  # deps
  require(ggplot2)
  require(dplyr)
  require(stringr)
  require(purrr)
  require(tibble)
  require(ggtree)     # Bioconductor
  require(patchwork)
  
  # Colors (match original)
  colors <- c(
    bp1_bp2 = "#1b9e77",
    bp1_bp3 = "#d95f02",
    bp2_bp3 = "#7570b3"
  )
  
  # 1) keys for this chromosome & comparisons
  keys <- grep(paste0("^", chrom, "__.*(bp1_bp2|bp1_bp3|bp2_bp3)"),
               names(twisst_data), value = TRUE)
  
  if (length(keys) == 0) stop("No matching keys for chromosome/comparisons.")
  
  # 2) summarize & pick one matching topology per comparison
  summary_list <- lapply(keys, function(k) {
    cmp <- sub(paste0(chrom, "__"), "", k)
    parts <- str_split(cmp, "_")[[1]]
    g1 <- toupper(parts[1]); g2 <- toupper(parts[2])
    tb <- twisst_data[[k]]
    
    keep <- which(sapply(tb$topologies, function(tr) {
      labs <- tr$tip.label
      (g1 %in% labs && g2 %in% labs) && ape::is.monophyletic(tr, c(g1, g2))
    }))
    
    if (length(keep) == 0) stop("No matching topology for ", cmp)
    
    wvec <- rowSums(tb$weights[, keep, drop = FALSE])
    list(
      comparison = cmp,
      data       = tibble(pos = tb$pos, weight = wvec),
      tree       = tb$topologies[[keep[1]]]
    )
  })
  
  # 3) combine & filter weights
  df_all <- bind_rows(lapply(summary_list, function(x) {
    cbind(comparison = x$comparison, x$data)
  })) %>%
    as_tibble() %>%
    filter(pos >= start, pos <= end) %>%
    mutate(comparison = factor(comparison, levels = names(colors)))
  
  # ---- plot builders ---------------------------------------------------------
  
  title_plot <- function(txt, col) {
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = txt, colour = col, size = 4) +
      theme_void() +
      theme(plot.margin = margin(0, 0, 0, 0))
  }
  
  tree_plot <- function(tr, col) {
    # ggtree respects plotting space; we keep it minimal and colored
    ggtree::ggtree(tr, size = 0.8, colour = col) +
      ggtree::geom_tiplab(colour = col, size = 3, hjust = 0, align = FALSE) +
      theme_void() +
      coord_cartesian(clip = "off") +
      theme(
        plot.margin = margin(2, 30, 2, 2),  # right margin for labels
      )
  }
  
  weights_plot <- function(df) {
    ggplot(df, aes(pos/1e6, weight, colour = comparison)) +
      geom_line(linewidth = 0.6) +
      geom_hline(yintercept = 1/3, linetype = "dashed", linewidth = 0.4, colour = "gray50") +
      scale_colour_manual(values = colors, breaks = names(colors),
                          labels = toupper(gsub("_", " vs ", names(colors)))) +
      coord_cartesian(ylim = c(0, 1)) +
      theme_classic(base_size = 11) +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0, face = "plain")
      )
  }
  
  # ---- assemble panels -------------------------------------------------------
  
  n <- length(summary_list)
  
  title_row <- map(summary_list, function(x) {
    title_text <- toupper(gsub("_", " vs ", x$comparison))
    title_plot(title_text, colors[[x$comparison]])
  })
  
  tree_row <- map(summary_list, function(x) {
    tree_plot(x$tree, colors[[x$comparison]])
  })
  
  bottom_plot <- weights_plot(df_all)+
    scale_x_continuous(expand = c(0.01, 0), name = "Chromosome Position (Mb)",  breaks = scales::breaks_extended(n = 4)) +
    theme_minimal() +
    theme(
      strip.background = element_blank(), 
      #strip.text = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.line.x = element_line(color = "black", size = 0.5),
      axis.ticks.x = element_line(color = "black", size = 0.5),
      axis.text.x = element_text(size = rel(2), family = "Arial", angle = 90, hjust = 1),
      axis.text.y = element_text(size = rel(2), family = "Arial", angle = 0, hjust = 1),
      axis.title.x = element_text(size = rel(2), family = "Arial", margin = margin(t = 20)),
      axis.title.y = element_text(size = rel(2), family = "Arial", margin = margin(r = 20)),
      legend.position = "right",
      legend.text = element_text(size=rel(1.2), family = "Arial"),
      panel.spacing = unit(2, "lines")
    ) +
    labs(y = expression("Topology weight"))
  
  # Layout: titles (row 1), trees (row 2), weights (row 3 spans all columns)
  # patchwork design string with n columns:
  top <- wrap_plots(title_row, ncol = n)
  mid <- wrap_plots(tree_row,  ncol = n)
  design <- "
  AAAAAAAAA
  BBBBBBBBB
  CCCCCCCCC
  "
  # repeat letters to exactly n columns in A/B, single C spanning n
  A <- wrap_plots(title_row, ncol = n)
  B <- wrap_plots(tree_row,  ncol = n)
  C <- bottom_plot
  
  # Use relative heights to mimic base plot layout (titles short, trees medium, weights tall)
  out <- (A / B / C) + plot_layout(heights = c(0.25 * tree_height, 1 * tree_height, 3))
  
  # Optionally scale overall width by tree_width (affects ggsave sizing more than on-screen)
  attr(out, "suggested_width")  <- max(6, n * 2.5 * tree_width)
  attr(out, "suggested_height") <- 0.25 * tree_height + 1 * tree_height + 3
  
  return(bottom_plot)
}

