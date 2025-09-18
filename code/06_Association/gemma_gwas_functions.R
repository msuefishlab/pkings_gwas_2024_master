library(ggrastr)

# functions inspired by erik enbody, with modifications

##USED##
plink.order <- function(plink.path, input.var){
  df.plink <- fread(plink.path)
  
  # Parse chr & position from SNP like "chr1_993362_T_C"
  snp_parts <- str_split_fixed(df.plink$SNP, "_", 4)
  
  df.plink <- df.plink %>%
    mutate(
      CHR_from_snp = snp_parts[, 1],
      ps           = as.numeric(snp_parts[, 2]),
      numeric_chr  = as.numeric(str_replace_all(CHR_from_snp, "chr|ups", ""))  # "chr1"→1
    ) %>%
    arrange(numeric_chr, ps) %>%
    mutate(
      CHR    = factor(CHR_from_snp),
      log_p  = -log10(.data[[input.var]]),
      chr_labels = CHR
    )
  
  # Row index shared across TESTs: 1,2,3,... in the order of unique (numeric_chr, ps)
  key <- paste(df.plink$numeric_chr, df.plink$ps, sep=":")
  df.plink$row <- as.integer(factor(key, levels = unique(key)))
  
  # return plain data.frame to match your other helpers
  df.plink <- as.data.frame(df.plink)
  df.plink
}

gemma.order <- function(gemma.path, input.var){
  df.gemma <- fread(gemma.path)
  df.gemma <- df.gemma %>% mutate(numeric_chr = as.numeric(str_replace_all(chr, "chr|ups", "")))
  df.gemma <- df.gemma %>% dplyr::arrange(numeric_chr, ps)
  df.gemma$chr <- factor(df.gemma$chr)
  df.gemma$row<-1:nrow(df.gemma)
  df.gemma <- as.data.frame(df.gemma)
  df.gemma$log_p <- -log10(df.gemma[,input.var])
  df.gemma$chr_labels <- df.gemma$chr
  chr_breaks <- df.gemma %>% group_by(chr) %>% dplyr::summarise(chr_breaks = mean(row))
  df.gemma
}

## USED##

peak_list_permutation_plink <- function(input.df, input.var, thresholdA, num_snps, test){
  autosome_peaks <- input.df %>% filter(TEST == test) %>% filter(!!sym(input.var) > thresholdA)
  
  #create granges object for individual SNPs that are above the threshold
  range.autosome_peaks <- GRanges(autosome_peaks$chr, IRanges(as.numeric(autosome_peaks$ps), as.numeric(autosome_peaks$ps)))
  
  #merge windows by some distance to other SNPs. defining it here by 75kb
  reduce.range.autosome_peaks <- GenomicRanges::reduce(range.autosome_peaks, min.gapwidth = 75000) #merge closest window
  
  #remove ranges that are made up of only less than num_snps 
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
  
  #give each "peak" a name
  names(reduce.range.autosome_peaks) <- 1:length(reduce.range.autosome_peaks)
  
  datalist <- list()
  
  #create a dataframe where each SNP gets an peak ID
  for (i in (1:length(reduce.range.autosome_peaks))){
    x<-autosome_peaks[queryHits(findOverlaps(range.autosome_peaks, reduce.range.autosome_peaks[i, ])), ]
    x$peak <- i
    datalist[[i]] <- x # add it to your list
  }
  peaks.named = do.call(rbind, datalist)
  
  return(peaks.named)
}

peak_list_permutation <- function(input.df, input.var, thresholdA, num_snps){
  autosome_peaks <- input.df %>% filter(!!sym(input.var) > thresholdA)
  
  #create granges object for individual SNPs that are above the threshold
  range.autosome_peaks <- GRanges(autosome_peaks$chr, IRanges(as.numeric(autosome_peaks$ps), as.numeric(autosome_peaks$ps)))
  
  #merge windows by some distance to other SNPs. defining it here by 75kb
  reduce.range.autosome_peaks <- GenomicRanges::reduce(range.autosome_peaks, min.gapwidth = 75000) #merge closest window
  
  #remove ranges that are made up of only less than num_snps 
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
  
  #give each "peak" a name
  names(reduce.range.autosome_peaks) <- 1:length(reduce.range.autosome_peaks)
  
  datalist <- list()
  
  #create a dataframe where each SNP gets an peak ID
  for (i in (1:length(reduce.range.autosome_peaks))){
    x<-autosome_peaks[queryHits(findOverlaps(range.autosome_peaks, reduce.range.autosome_peaks[i, ])), ]
    x$peak <- i
    datalist[[i]] <- x # add it to your list
  }
  peaks.named = do.call(rbind, datalist)
  
  return(peaks.named)
}


plink.peak.bed <- function(peaks.named, log_p_threshold){
  peaks.named %>%
    mutate(
      ps    = as.numeric(ps),
      log_p = as.numeric(log_p)
    ) %>%
    group_by(numeric_chr, peak) %>%
    summarise(
      start_bp = min(ps, na.rm = TRUE),
      end_bp   = max(ps, na.rm = TRUE),
      region_len = max(ps, na.rm = TRUE) - min(ps, na.rm = TRUE),  # avoid start/end refs
      num_snps  = n(),
      mean_log_p = mean(log_p, na.rm = TRUE),
      max_log_p  = max(log_p, na.rm = TRUE),
      num_snps_above_threshold = sum(log_p > log_p_threshold, na.rm = TRUE),
      most_significant_bp = ps[which.max(log_p)],
      .groups = "drop"
    )
}

###USED###
gemma.peak.bed <- function(peaks.named, log_p_threshold){
  # Create a bed file that has a summary of each of these regions
  peak.bed <- peaks.named %>% 
    group_by(numeric_chr, peak) %>% 
    summarise(start = min(ps),
              end = max(ps),
              length = end - start,
              num.snps = n(),
              mean.log_p = mean(log_p),
              max.log_p = max(log_p),
              num.snps.above.threshold = sum(log_p > log_p_threshold),
              most_significant_bp = ps[which.max(log_p)],  # Get the position of the most significant SNP
              .groups = "keep")
  return(peak.bed)
}

### USED###
gemma.peak.ROW.bed <- function(peaks.named){
  #create a bed file that has a summary of each of these regions
  peak.bed <- peaks.named %>% group_by(chr,peak) %>%
    summarise(start = min(row),
              end = max(row),
              length = end - start,
              num.snps = n(),
              mean.log_p = mean(log_p),
              max.log_p = max(log_p),
              .groups = "keep")
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
    chr      = unique(df.in$chr),
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
    chr      = unique(df.in$chr),
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
