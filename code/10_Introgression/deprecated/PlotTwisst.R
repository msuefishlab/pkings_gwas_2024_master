#!/usr/bin/env Rscript

#Load twisst plotting functions
source("twisst/plot_twisst.R")

#Load required libraries
library(dplyr)
library(ggplot2)
library(ggthemes)

#weights file with a column for each topology
weights_file <- c("chr1_100snps/chr1.phyml_bionj.w100snps.weights.csv.gz",
                  "chr1A_100snps/chr1A.phyml_bionj.w100snps.weights.csv.gz",
                  "chr2_100snps/chr2.phyml_bionj.w100snps.weights.csv.gz",
                  "chr3_100snps/chr3.phyml_bionj.w100snps.weights.csv.gz",
                  "chr5_100snps/chr5.phyml_bionj.w100snps.weights.csv.gz",
                  "chr7_100snps/chr7.phyml_bionj.w100snps.weights.csv.gz",
                  "chr9_100snps/chr9.phyml_bionj.w100snps.weights.csv.gz",
                  "chr25_100snps/chr25.phyml_bionj.w100snps.weights.csv.gz")

#coordinates file for each window
window_data_file <- c("chr1_100snps/chr1.phyml_bionj.w100snps.data.tsv",
                      "chr1A_100snps/chr1A.phyml_bionj.w100snps.data.tsv",
                      "chr2_100snps/chr2.phyml_bionj.w100snps.data.tsv",
                      "chr3_100snps/chr3.phyml_bionj.w100snps.data.tsv",
                      "chr5_100snps/chr5.phyml_bionj.w100snps.data.tsv",
                      "chr7_100snps/chr7.phyml_bionj.w100snps.data.tsv",
                      "chr9_100snps/chr9.phyml_bionj.w100snps.data.tsv",
                      "chr25_100snps/chr25.phyml_bionj.w100snps.data.tsv")

#Load data using import.twisst function
twisst_data <- import.twisst(weights_files=weights_file, window_data_files=window_data_file)

#Reorder topos according to mean weight 
topo_order <- order(twisst_data$weights_overall_mean, decreasing=T)
twisst_data <- subset.twisst.by.topos(twisst_data, topo_order)

#Set P.inornata as root
for (i in 1:length(twisst_data$topos)) twisst_data$topos[[i]] <- root(twisst_data$topos[[i]], "P.inornata", resolve.root = T)

#Extract weights for each chromosome, and append zeros to start/end to allow
#polygon drawing (in the end we didnt use polygons for plotting but this remains here
#as a legacy to an idea)
chrom.list <- c("chr1","chr1A","chr2","chr3","chr5","chr7","chr9","chr25")
for (CHROM in 1:length(chrom.list)){
  temp <- data.frame(twisst_data$pos[CHROM], twisst_data$weights[CHROM])
  colnames(temp) <- c("pos", "topo1", "topo2", "topo3")
  temp$chr <- chrom.list[CHROM]
  top <- head(temp,1)
  top$pos <- top$pos-1
  top$topo1 <- 0
  top$topo2 <- 0
  top$topo3 <- 0
  tail <- tail(temp,1)
  tail$pos <- tail$pos-1
  tail$topo1 <- 0
  tail$topo2 <- 0
  tail$topo3 <- 0
  temp <- rbind(top, temp, tail)
  assign(paste0("weights_",chrom.list[CHROM]), temp)
  rm(top, temp, tail)
}

#combine into  single data.frame
weights_all_chroms <- rbind(weights_chr1,weights_chr1A,weights_chr2,weights_chr3,
                            weights_chr5,weights_chr7,weights_chr9,weights_chr25)

#Load 28 loci positions and calculate centre
Loci28Pos <- read.delim("28Loci_start_stop.txt")
Loci28Pos$mid <- (Loci28Pos$start.pos + Loci28Pos$end.pos)/2
Loci28Pos$chr <- as.factor(Loci28Pos$chr)

#Plot topologies
topo_cols <- c("#FF0010", #Red
               "#00998F", #Turquoise
               "#003380") #Navy
pdf("plots/TWISST_100snps_trees.pdf", width=5, height=2)
par(mfrow = c(1,length(twisst_data$topos)), mar = c(1,1,2,1), xpd=NA)
for (n in 1:length(twisst_data$topos)){
  plot.phylo(twisst_data$topos[[n]], type = "unrooted", edge.color=topo_cols[n], edge.width=5, label.offset=.1, cex = 1)
  mtext(side=3,text=paste0("topo",n))
}
dev.off()

#Plot weighting across chroms
#Reorder factor levels so that chroms are plotted in the desired order
weights_all_chroms$chr <- as.factor(weights_all_chroms$chr)
weights_all_chroms$chr <- factor(weights_all_chroms$chr, levels=c("chr1","chr1A","chr2","chr3","chr5","chr7","chr9","chr25"))

#Plot them ...
pdf("plots/TWISST_100snps_all_28locichroms.pdf", height = 14, width = 16)
print(
  ggplot(data = weights_all_chroms) +
    geom_rect(data = Loci28Pos, aes(xmin = start.pos/1000000, xmax = end.pos/1000000, ymin = 0, ymax = 1.1), fill = "#D3D3D3") +
    geom_text(data = Loci28Pos, aes(x = mid/1000000, y = 1.2, label = loci.no), size = 2.5) +
    geom_col(aes(x=pos/1000000, y=topo1), col = topo_cols[1], size = 0.25) +
    geom_col(aes(x=pos/1000000, y=topo2), col = topo_cols[2], size = 0.25) +
    geom_col(aes(x=pos/1000000, y=topo3), col = topo_cols[3], size = 0.25) +
    scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0)) +
    scale_x_continuous(limits = c(-2,max(weights_all_chroms$pos/1000000)), expand = c(0, 0)) +
    ylab("Weighting") +
    xlab("Mb") +
    facet_grid(chr~., scales = "free_x", space = "free_x"))
dev.off()


#Create violin / box plots showing topology 1 weighting
#Determine which windows are genomic background and which are wihtin
#the 28 loci
weights_all_chroms$status <- "Genomic background"
for (LOCI_NO in 1:28){
  for (ROW in 1:nrow(weights_all_chroms)){
    if (weights_all_chroms$chr[ROW] == Loci28Pos$chr[LOCI_NO] & 
        between(weights_all_chroms$pos[ROW], Loci28Pos$start.pos[LOCI_NO], Loci28Pos$end.pos[LOCI_NO]) == TRUE){
      weights_all_chroms$status[ROW] <- "28 loci"
    }
  } 
}
weights_all_chroms$status <- as.factor(weights_all_chroms$status)
weights_all_chroms$status <- factor(weights_all_chroms$status, levels=c("Genomic background","28 loci"))

#Plot ...
pdf("plots/TWISST_100snps_boxplots_topo1.pdf", width=4.5, height=4)
library(ggpubr)
print(ggviolin(weights_all_chroms, x = "status", y = "topo1", fill = "status",
         palette = c("#CACACA", "#CACACA"),
         add = "boxplot", add.params = list(fill = "white")) +
  ylab("Topology 1 weight") +
  xlab("") +
  theme(legend.position = "none"))
dev.off()


#Create zoom-in plot of chromosome 1A
chr1A_Loci <- Loci28Pos %>% filter(chr == "chr1A")
pdf("plots/TWISST_100snps_chr1A.pdf", height = 3, width = 12)
ggplot(data = weights_chr1A) +
  geom_rect(data = chr1A_Loci, aes(xmin = start.pos/1000000, xmax = end.pos/1000000, ymin = 0, ymax = 1.1), fill = "#D3D3D3") +
  geom_text(data = chr1A_Loci, aes(x = mid/1000000, y = 1.2, label = loci.no), size = 3) +
  geom_col(aes(x=pos/1000000, y=topo1), col = topo_cols[1], size = 0.5) +
  geom_col(aes(x=pos/1000000, y=topo2), col = topo_cols[2], size = 0.5) +
  geom_col(aes(x=pos/1000000, y=topo3), col = topo_cols[3], size = 0.5) +
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0)) +
  scale_x_continuous(limits = c(-1,max(weights_chr1A$pos/1000000)), expand = c(0, 0)) +
  ylab("Weighting") +
  xlab("Chr1A (Mb)")
dev.off()