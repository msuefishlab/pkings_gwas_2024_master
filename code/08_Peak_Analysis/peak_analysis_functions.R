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

get_stat_values<-function(df, start, end, population) {
  df %>%
    filter(
      window_pos_2 >= start,
      window_pos_1 <= end,
      pop1 == population
    ) %>%
    pull(value)
}

get_data4plot <- function(df.bed, df, stat) {
  # Initialize a list to store results
  result_list <- list()
  
  # Get unique populations from the dataset
  populations <- unique(df$pop1)
  
  # Loop through each peak
  for (i in 1:nrow(df.bed)) {
    peak <- df.bed[i, ]
    chr <- paste0("chr", peak$chr)
    start <- peak$start
    end <- peak$end
    lg <- peak$lg
    
    # Pre-filter df for the current chromosome
    df_chr <- df %>% filter(chromosome == chr)
    
    ## Process the peak region ##
    for (pop in populations) {
      # Obtain the statistic values for the peak region
      peak_stat <- get_stat_values(df_chr, start, end, pop)
      
      # Only add to result if there is at least one value
      if (length(peak_stat) > 0) {
        result_list[[length(result_list) + 1]] <- data.frame(
          peak = rep(i, length(peak_stat)),        # Peak number repeated
          population = rep(pop, length(peak_stat)),  # Population name repeated
          interval = rep(paste("peak", i, sep = "_"), length(peak_stat)),
          stringsAsFactors = FALSE
        ) %>% mutate(!!stat := peak_stat)  # Dynamically name the statistic column
      }
    }
    
    ## Process random intervals ##
    # Determine the maximum possible position for random interval selection
    max_pos <- max(df_chr$window_pos_2)
    for (pop in populations) {
      # Generate a data frame of random intervals (100 per population for the current peak)
      random_intervals <- data.frame(
        thestart = sample(1:(max_pos - lg), 100)
      ) %>%
        mutate(
          theend = thestart + lg  # Create the end position of each interval
        ) %>%
        filter(theend <= max_pos, thestart < theend)  # Exclude invalid intervals
      
      # Initialize a list to store data frames for each random interval
      random_dfs <- list()
      for (j in 1:nrow(random_intervals)) {
        # For each random interval, get the statistic values
        stat_values <- get_stat_values(df_chr, random_intervals$thestart[j], random_intervals$theend[j], pop)
        # Only create a data frame if stat_values is non-empty
        if (length(stat_values) > 0) {
          random_dfs[[length(random_dfs) + 1]] <- data.frame(
            peak = as.factor(i),
            population = rep(paste0(pop, "X"), length(stat_values)),
            interval = rep(paste("random", i, j, sep = "_"), length(stat_values)),
            stringsAsFactors = FALSE
          ) %>% mutate(!!stat := stat_values)
        }
      }
      
      # Only add the combined random interval data if there is any
      if (length(random_dfs) > 0) {
        result_list[[length(result_list) + 1]] <- do.call(rbind, random_dfs)
      }
    }
  }
  
  # Combine all results into a single data frame
  result_data <- do.call(rbind, result_list)
  
  # Dynamically filter out any NA values based on the stat column
  result_data <- result_data %>%
    filter(!is.na(!!sym(stat)))
  
  return(result_data)
}



get_fst_stat_values<-function(df, start, end, comparison) {
  df %>%
    filter(
      window_pos_2 >= start,
      window_pos_1 <= end,
      comparison == !!comparison
    ) %>%
    pull(value)
}

get_fst_stat_values <- function(df, start, end, comparison) {
  df %>%
    filter(
      window_pos_2 >= start,
      window_pos_1 <= end,
      comparison == !!comparison
    ) %>%
    pull(value)
}

get_fst_data4plot <- function(df.bed, df, stat = "avg_wc_fst") {
  # Initialize a list to store results
  result_list <- list()
  
  # Get unique comparisons from the dataset
  comparisons <- unique(df$comparison)
  
  # Loop through each peak
  for (i in 1:nrow(df.bed)) {
    peak <- df.bed[i, ]
    chr <- paste0("chr", peak$chr)
    start <- peak$start
    end <- peak$end
    lg <- peak$lg
    
    # Pre-filter df for the current chromosome
    df_chr <- df %>% filter(chromosome == chr)
    
    ## Process the peak region ##
    for (comp in comparisons) {
      # Obtain the FST values for the peak region
      peak_stat <- get_fst_stat_values(df_chr, start, end, comp)
      
      # Only add to results if any values were found
      if (length(peak_stat) > 0) {
        result_list[[length(result_list) + 1]] <- data.frame(
          peak = rep(i, length(peak_stat)),
          comparison = rep(comp, length(peak_stat)),
          interval = rep(paste("peak", i, sep = "_"), length(peak_stat)),
          stringsAsFactors = FALSE
        ) %>% mutate(!!stat := peak_stat)
      }
    }
    
    ## Process random intervals ##
    # Determine the maximum possible position for random interval selection
    max_pos <- max(df_chr$window_pos_2)
    for (comp in comparisons) {
      # Generate a data frame of 100 random intervals for the current peak and comparison
      random_intervals <- data.frame(
        thestart = sample(1:(max_pos - lg), 100)
      ) %>%
        mutate(
          theend = thestart + lg
        ) %>%
        filter(theend <= max_pos, thestart < theend)
      
      # Loop over each random interval so we can tag its origin
      for (j in 1:nrow(random_intervals)) {
        stat_values <- get_fst_stat_values(
          df_chr, 
          random_intervals$thestart[j], 
          random_intervals$theend[j], 
          comp
        )
        # Only add to results if there are statistic values
        if (length(stat_values) > 0) {
          result_list[[length(result_list) + 1]] <- data.frame(
            peak = rep(i, length(stat_values)),
            comparison = rep(paste0(comp, "X"), length(stat_values)),
            interval = rep(paste("random", i, j, sep = "_"), length(stat_values)),
            stringsAsFactors = FALSE
          ) %>% mutate(!!stat := stat_values)
        }
      }
    }
  }
  
  # Combine all results into a single data frame
  result_data <- do.call(rbind, result_list)
  
  # Dynamically filter out any NA values based on the stat column
  result_data <- result_data %>%
    filter(!is.na(!!sym(stat)))
  
  return(result_data)
}

generate_palette_alpha <- function(pops, phenos, phenotype_palettes) {
  # Initialize an empty palette
  palette <- c()
  
  # Get unique phenotypes and their associated populations
  unique_phenos <- unique(phenos)
  
  for (pheno in unique_phenos) {
    # Get the Brewer palette for the phenotype
    brewer_palette <- phenotype_palettes[[pheno]]
    if (is.null(brewer_palette)) {
      stop(paste("Phenotype", pheno, "not recognized!"))
    }
    
    # Get populations belonging to this phenotype
    group_indices <- which(phenos == pheno)
    group_pops <- pops[group_indices]
    num_pops <- length(group_pops)
    
    # Generate a palette with enough distinct shades, reversed
    base_colors <- rev(brewer.pal(min(max(3, num_pops), 9), brewer_palette)) # Reverse for darkest first
    
    # If more populations than colors, interpolate additional colors
    if (num_pops > length(base_colors)) {
      base_colors <- rev(colorRampPalette(base_colors)(num_pops)) # Reverse interpolated palette
    }
    
    # Assign colors for each population
    for (i in seq_along(group_pops)) {
      pop <- group_pops[i]
      
      # Dark color for the main population
      color <- base_colors[i]
      
      # Assign the same color for "R" but with reduced alpha
      palette[pop] <- color
      palette[paste0(pop, "R")] <- scales::alpha(color, 0.4) # Lower opacity
    }
  }
  
  return(palette)
}


# Function to generate FST palette with alpha
generate_fst_palette_alpha <- function(comparisons, compare_colors) {
  # Initialize an empty palette
  palette <- c()
  
  # Get unique compare_color groups
  unique_compare_colors <- unique(comparisons$compare_color)
  
  for (compare_color in unique_compare_colors) {
    # Get the Brewer palette for the compare_color group
    brewer_palette <- compare_colors[[compare_color]]
    if (is.null(brewer_palette)) {
      stop(paste("Compare_color", compare_color, "not recognized!"))
    }
    
    # Get rows corresponding to this compare_color
    rows <- which(comparisons$compare_color == compare_color)
    num_rows <- length(rows)
    
    # Generate a palette with enough distinct shades, reversed
    base_colors <- rev(brewer.pal(min(max(3, num_rows), 9), brewer_palette)) # Reverse for darkest first
    
    # If more comparisons than colors, interpolate additional colors
    if (num_rows > length(base_colors)) {
      base_colors <- rev(colorRampPalette(base_colors)(num_rows)) # Reverse interpolated palette
    }
    
    # Assign colors for each comparison
    for (i in seq_along(rows)) {
      comparison <- comparisons$comparison[rows[i]]
      
      # Dark color for the main comparison
      color <- base_colors[i]
      
      # Assign the same color for "R" but with reduced alpha
      palette[comparison] <- color
      palette[paste0(comparison, "R")] <- scales::alpha(color, 0.4)
    }
  }
  
  return(palette)
}