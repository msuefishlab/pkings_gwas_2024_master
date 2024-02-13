# functions inspired by erik enbody, with modifications

gemma.order <- function(gemma.path, input.var){
  df.gemma <- fread(gemma.path)
  df.gemma <- df.gemma %>% mutate(numeric_chr = as.numeric(str_replace_all(chr, "ups", "")))
  df.gemma <- df.gemma %>% dplyr::arrange(numeric_chr, ps)
  df.gemma$chr <- factor(df.gemma$chr)
  df.gemma$row<-1:nrow(df.gemma)
  df.gemma <- as.data.frame(df.gemma)
  df.gemma$log_p <- -log10(df.gemma[,input.var])
  df.gemma$df.log_p_rollmean <- zoo::rollmean(df.gemma[,input.var],50,fill=NA)
  df.gemma$chr_labels <- df.gemma$chr
  chr_breaks <- df.gemma %>% group_by(chr) %>% dplyr::summarise(chr_breaks = mean(row))
  df.gemma
}

peak_list_permutation <- function(input.df, input.var, thresholdA, thresholdZ,num_snps){
  input.df <- as.data.frame(input.df)
  
  #                                                          %>% 
  #                                                            mutate(chr_ordered = factor(chr_ordered, levels = chr_order)) %>%
  #                                                            arrange(chr_ordered)
  
  auto_fst <- input.df %>% filter(!grepl("ups", chr))

  peaks_fst <- auto_fst %>% filter(!!sym(input.var) > thresholdA)
  
  #fst_range <- GRanges(peaks_fst$chr, IRanges(as.numeric(peaks_fst$ps), as.numeric(peaks_fst$ps)))
  
  #create granges object for individual SNPs that are above the threshold
  range.peaks_fst <- GRanges(peaks_fst$chr, IRanges(as.numeric(peaks_fst$ps), as.numeric(peaks_fst$ps)))
  
  #merge windows by some distance to other SNPs. defining it here by 75kb
  reduce.range.peaks_fst <- GenomicRanges::reduce(range.peaks_fst, min.gapwidth = 75000) #merge closest window
  #remove ranges that are made up of only 1 SNP. 
  reduce.range.peaks_fst <- reduce.range.peaks_fst[width(reduce.range.peaks_fst) > num_snps]
  #give each "peak" a name
  names(reduce.range.peaks_fst) <- 1:length(reduce.range.peaks_fst)
  
  datalist <- list()
  #create a dataframe where each SNP gets an peak ID
  for (i in (1:length(reduce.range.peaks_fst))){
    x<-peaks_fst[queryHits(findOverlaps(range.peaks_fst, reduce.range.peaks_fst[i, ])), ]
    x$peak <- i
    datalist[[i]] <- x # add it to your list
  }
  peaks.named = do.call(rbind, datalist)
  
  #an alternative way to make the bed file directly from granges, but it comes with 0 metadata
  #peak.bed <- as.data.frame(reduce.range.peaks_fst)
  #peak.bed$peak <- names(reduce.range.peaks_fst)
  
  return(peaks.named)
  
}

gemma.peak.bed <- function(peaks.named){
  #create a bed file that has a summary of each of these regions
  peak.bed <- peaks.named %>% group_by(chr,peak) %>% 
    summarise(start = min(ps),
              end = max(ps),
              length = end - start,
              num.snps = n(),
              mean.log_p = mean(log_p),
              max.log_p = max(log_p),
              .groups = "keep")
  return(peak.bed)
}

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

gemma.get_genes <- function(comp.bed){
  parv.df.an.gr <- GRanges(parv.df.an$seqid, IRanges(as.numeric(parv.df.an$start), as.numeric(parv.df.an$end)))
  
  comp.bed <- as.data.frame(comp.bed)
  comp.bed.gr <- GRanges(comp.bed$chr, IRanges(comp.bed$start, comp.bed$end), peak = comp.bed$peak)
  
  comp_overlap <- findOverlaps(comp.bed.gr, parv.df.an.gr)
  
  comp.annotated <- cbind(comp.bed[queryHits(comp_overlap), ], parv.df.an[subjectHits(comp_overlap), c(-1, -2, -3)])
  
  comp.annotated <- comp.annotated %>% mutate(gene_name = gene_symbol) %>% 
    select(-gene_symbol, -Name)
  
  comp.annotated
}

manc2 <- function(df.in, input.var){
  chr_breaks <- df.in %>% filter(chr_ordered != "chrunknown" & !is.na(row)) %>% 
    mutate(chr_ordered = factor(chr_ordered)) %>%
    group_by(chr_ordered, chr_labels) %>% 
    dplyr::summarise(chr_breaks = mean(row))
  
  chrom.colors <- data.frame(chr_ordered=unique(df.in$chr_ordered),
                             color.num = rep(1:2,length(unique(df.in$chr_ordered)))) %>% 
    distinct(chr_ordered, .keep_all = T)
  
  df.in2 <- df.in %>% #mutate(row = 1:n()) %>% 
    left_join(chrom.colors, by = "chr_ordered") %>% 
    mutate(color.num = as.factor(color.num))
  
  df.in2 %>% filter(chr_labels != "unknown" & !is.na(row)) %>%
    ggplot(aes_string(x = "row", y = input.var, col = "color.num")) + theme_bw() +
    theme(legend.position="none",
          #panel.border=element_blank(),
          panel.border = element_blank(), axis.line.x = element_line(), axis.line.y = element_line(),
          axis.title.x=element_blank(),
          #axis.text.x = element_text(angle = 45, color = "black"),
          axis.text.x = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor.y=element_blank(),
          axis.title.y = element_text(size=6),
          axis.text = element_text(size=6),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(size=0.2)) +
    geom_point(size=0.9,shape=20,stroke=0.2) +
    scale_color_manual(values=rep(c("grey30","grey70"))) +
    #scale_color_manual(values = c(rep_len(c("grey30", "red"), length(unique(chr_breaks$chr_ordered))+1))) #+
    #scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1),minor_breaks = NULL) +
    scale_x_continuous(breaks = chr_breaks$chr_breaks, 
                       labels = function(labels) {
                         sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n'), chr_breaks$chr_labels[i]))
                       }) +
    labs(y=expression(paste("-log"[10], italic("P"),"-value"))) 
  #scale_x_continuous(breaks=chr_breaks$chr_breaks, 
  #                   labels = chr_breaks$chr_labels)
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