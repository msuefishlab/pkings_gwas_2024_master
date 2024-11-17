# Define functions for each common task

# Function to filter data
filter_data <- function(data, chr, start, end, window_size = 200000) {
  data %>% filter(chromosome == chr & window_pos_1 >= start - window_size & window_pos_1 <= end + window_size)
}

# Function to create GWAS plot
create_gwas_plot <- function(pk.lmm, filtered_gwas_ld_data, idx_snp_peak, overlapping_annotations_df, chr.idx, start,end,window_size=200000) {
  
  # Filter pk.lmm to get p_lrt values within the specified range
  min_p_lrt <- min(pk.lmm$p_lrt[pk.lmm$ps >= start & pk.lmm$ps <= end], na.rm = TRUE)
  
  # Calculate the corresponding -log10(min_p_lrt) for y-axis limit
  y_min <- -log10(min_p_lrt)
  
  ggplot() +
    geom_point(data = pk.lmm, mapping = aes(x = ps, y = -log10(p_lrt)), color = "grey") +
    geom_point(data = filtered_gwas_ld_data, mapping = aes(x = BP_B, y = -log10(p_lrt), color = R2)) +
    scale_color_gradientn(colours = c("blue", "green", "yellow", "red")) +
    geom_hline(yintercept = significant, linetype = "dashed", color = "red", size = 0.5) +
    geom_hline(yintercept = suggestive, linetype = "dashed", color = "blue", size = 0.5) +
    geom_point(data = idx_snp_peak, mapping = aes(x = BP_A, y = -log10(p_lrt)), shape = 5, size = 5, fill = "black") +
    geom_gene_arrow(data = overlapping_annotations_df, mapping = aes(xmin = start, xmax = end, y = -1, fill = external_gene_name), show.legend = FALSE) +
    geom_gene_label(data = overlapping_annotations_df, mapping = aes(xmin = start, xmax = end, y = -1, label = external_gene_name), align = "left", show.legend = FALSE) +
    xlim(start - window_size, end + window_size) +
    ylim(NA, y_min) +  # Set the y-axis lower limit dynamically; upper limit will adjust automatically
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(size = rel(1.2), family = "Arial", angle = 25, hjust = 1),
      axis.text.y = element_text(size = rel(1.2), family = "Arial", angle = 0, hjust = 1),
      axis.title.x = element_text(size = rel(1.2), family = "Arial", margin = margin(t = 20)),
      axis.title.y = element_text(size = rel(1.2), family = "Arial", margin = margin(r = 20)),
      legend.position = "right"
    ) +
    xlab(paste0("position along ", chr.idx)) +
    ylab("-log(P)")
}

# Function to create nucleotide diversity (π) plot
create_pi_plot <- function(pi_peak, df.bed.idx, chr.idx) {
  # Generate a color palette
  n_bp_populations <- dim(as.data.frame(pi_peak %>% filter(phenotype_pop1=="BP") %>% ungroup()%>% dplyr::select(pop1) %>% distinct()))[1]
  n_tp_populations <-dim(as.data.frame(pi_peak %>% filter(phenotype_pop1=="TP") %>% ungroup()%>% dplyr::select(pop1) %>% distinct()))[1]
  bp_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "Blues")[4:9])(n_bp_populations)  # Darker blues
  tp_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "Reds")[4:9])(n_tp_populations)  # Darker reds
  
  # Create a named vector for manual color mapping
  color_mapping <- c(
    setNames(bp_colors, unique(pi_peak$pop1[pi_peak$phenotype_pop1 == "BP"])),
    setNames(tp_colors, unique(pi_peak$pop1[pi_peak$phenotype_pop1 == "TP"]))
  )
  
  # Add a combined key to the data
  pi_peak$combined_key <- paste0(pi_peak$phenotype_pop1, "_", pi_peak$pop1)
  
  ggplot() +
    geom_line(
      data = pi_peak,
      aes(
        x = window_pos_1,
        y = value,
        colour = pop1  # Map the combined key to the color
      )
    ) +
    geom_rect(
      data = df.bed.idx,
      aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
      fill = "gray", alpha = 0.2
    ) +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(size = rel(1.2), family = "Arial", angle = 25, hjust = 1),
      axis.text.y = element_text(size = rel(1.2), family = "Arial", angle = 0, hjust = 1),
      axis.title.x = element_text(size = rel(1.2), family = "Arial", margin = margin(t = 20)),
      axis.title.y = element_text(size = rel(1.2), family = "Arial", margin = margin(r = 20)),
      legend.position = "right"
    ) +
    scale_y_continuous(limits = c(0, .01), labels = percent_format()) +
    scale_colour_manual(
      values = color_mapping  # Apply the color mapping
    ) +
    xlab(paste0("position along ", chr.idx)) +
    ylab("π (% Nucleotide Diversity)") +
    facet_wrap(~ phenotype_pop1,ncol = 1)
}

# Function to create Tajima's D plot
create_tajima_plot <- function(tpp_plot, df.bed.idx, chr.idx) {
  
  # Generate a color palette
  n_bp_populations <- dim(as.data.frame(tpp_plot %>% filter(phenotype_pop1=="BP") %>% ungroup()%>% dplyr::select(pop1) %>% distinct()))[1]
  n_tp_populations <-dim(as.data.frame(tpp_plot %>% filter(phenotype_pop1=="TP") %>% ungroup()%>% dplyr::select(pop1) %>% distinct()))[1]
  bp_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "Blues")[4:9])(n_bp_populations)  # Darker blues
  tp_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "Reds")[4:9])(n_tp_populations)  # Darker reds
  
  # Create a named vector for manual color mapping
  color_mapping <- c(
    setNames(bp_colors, unique(tpp_plot$pop1[tpp_plot$phenotype_pop1 == "BP"])),
    setNames(tp_colors, unique(tpp_plot$pop1[tpp_plot$phenotype_pop1 == "TP"]))
  )
  
  # Add a combined key to the data
  tpp_plot$combined_key <- paste0(tpp_plot$phenotype_pop1, "_", tpp_plot$pop1)
  
  ggplot() +
    geom_line(data = tpp_plot, aes(x = window_pos_1, y = value, colour = pop1)) +
    geom_rect(
      data = df.bed.idx,
      aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
      fill = "gray", alpha = 0.2
    ) +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(size = rel(1.2), family = "Arial", angle = 25, hjust = 1),
      axis.text.y = element_text(size = rel(1.2), family = "Arial", angle = 0, hjust = 1),
      axis.title.x = element_text(size = rel(1.2), family = "Arial", margin = margin(t = 20)),
      axis.title.y = element_text(size = rel(1.2), family = "Arial", margin = margin(r = 20)),
      legend.position = "right"
    ) +
    scale_y_continuous(limits = c(-2, 4), breaks = c(-2, 0, 2, 4)) +
    scale_colour_manual(
      values = color_mapping  # Apply the color mapping
    ) +
    xlab(paste0("position along ", chr.idx)) +
    ylab("Tajima's D") +
    facet_wrap(~ phenotype_pop1,ncol = 1)
}

# Function to create Fst plot
create_fst_plot <- function(fpp_plot, df.bed.idx, chr.idx) {
  
  # Generate a color palette
  same_bp_populations <- dim(as.data.frame(fpp_plot %>% filter(compare_color=="BP-BP") %>% ungroup()%>% dplyr::select(comparison) %>% distinct()))[1]
  same_tp_populations <-dim(as.data.frame(fpp_plot %>% filter(compare_color=="TP-TP") %>% ungroup()%>% dplyr::select(comparison) %>% distinct()))[1]
  different_populations <- dim(as.data.frame(fpp_plot %>% filter(compare_color=="BP-TP") %>% ungroup()%>% dplyr::select(comparison) %>% distinct()))[1]
  
  bp_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "Blues")[4:9])(same_bp_populations)  # Darker blues
  tp_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "Reds")[4:9])(same_tp_populations)  # Darker reds
  mixed_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "Purples")[4:9])(different_populations)  # Darker purples
  
  
  # Create a named vector for manual color mapping
  color_mapping <- c(
    setNames(bp_colors, unique(fpp_plot$comparison[fpp_plot$compare_color == "BP-BP"])),
    setNames(tp_colors, unique(fpp_plot$comparison[fpp_plot$compare_color == "TP-TP"])),
    setNames(mixed_colors, unique(fpp_plot$comparison[fpp_plot$compare_color == "BP-TP"]))
  )
  
  ggplot() +
    geom_line(data = fpp_plot, aes(x = window_pos_1, y = value, colour = comparison)) +
    geom_rect(
      data = df.bed.idx,
      aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
      fill = "gray", alpha = 0.2
    ) +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text = element_text(size = rel(1.2), family = "Arial", angle = 25, hjust = 1),
      axis.title.x = element_text(size = rel(1.2), family = "Arial", margin = margin(t = 20)),
      axis.title.y = element_text(size = rel(1.2), family = "Arial", margin = margin(r = 20)),
      legend.position = "right"
    ) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    scale_colour_manual(
      values = color_mapping  # Apply the color mapping
    ) +
    ylim(0, 1) +
    xlab(paste0("position along ", chr.idx)) +
    ylab("Weir & Cockerham's Fst") +
    facet_wrap(~ comparison_classification, ncol=1)
}


