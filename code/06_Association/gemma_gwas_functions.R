# functions inspired by erik enbody, with modifications

##USED##
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
peak_list_permutation <- function(input.df, input.var, thresholdA, num_snps){
  input.df <- as.data.frame(input.df)
  autosome_peaks <- autosome_data %>% filter(!!sym(input.var) > thresholdA)
  
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

###USED###
gemma.peak.bed <- function(peaks.named){
  #create a bed file that has a summary of each of these regions
  peak.bed <- peaks.named %>% group_by(numeric_chr,peak) %>% 
    summarise(start = min(ps),
              end = max(ps),
              length = end - start,
              num.snps = n(),
              mean.log_p = mean(log_p),
              max.log_p = max(log_p),
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

manc2.labels <- function(df.in, input.var){
  # Assuming df.in is your input dataframe
  chr_breaks <- df.in %>%
    group_by(chr) %>% 
    dplyr::summarise(chr_breaks = mean(row))
  
  chrom.colors <- data.frame(chr=unique(df.in$chr),
                             color.num = rep(1:2,length(unique(df.in$chr)))) %>% 
    distinct(chr, .keep_all = T)
  
  df.in2 <- df.in %>%
    left_join(chrom.colors, by = "chr") %>% 
    mutate(color.num = as.factor(color.num))
  
  df.in2 %>%
    ggplot(aes_string(x = "row", y = input.var, col = "color.num",size=input.var)) + theme_bw() +
    theme(legend.position="none",
          #panel.border=element_blank(),
          panel.border = element_blank(), axis.line.x = element_line(), axis.line.y = element_line(),
          #axis.title.x=element_blank(),
          #axis.text.x = element_text(angle = 45, color = "black"),
          #axis.text.x = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor.y=element_blank(),
          axis.title.y = element_text(size=6),
          axis.text = element_text(size=6),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(size=0.2)) +
    geom_point(shape=20,stroke=0.2,alpha=0.25) +
    scale_color_manual(values=rep(c("grey30","grey70"))) +
    #scale_color_manual(values = c(rep_len(c("grey30", "red"), length(unique(chr_breaks$chr_ordered))+1))) #+
    #scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1),minor_breaks = NULL) +
    scale_x_continuous(breaks = chr_breaks$chr_breaks, 
                       labels = function(labels) {
                         sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n'), chr_breaks$chr[i]))
                       }) +
    labs(y=expression(paste("-log"[10], italic("P"),"-value"))) 
  #scale_x_continuous(breaks=chr_breaks$chr_breaks, 
  #                   labels = chr_breaks$chr_labels)
}