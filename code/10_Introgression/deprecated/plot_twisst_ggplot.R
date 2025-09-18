library(dplyr)
library(tidyr)
library(ggplot2)


#' Refactored ggplot2 version of plot.twisst
#' 
#' @param twisst_object List with $weights (list of data.frames), $pos or $interval_data, $topos
#' @param xlim 
plot.twisst.gg <- function(twisst_object,chr,
                        xlim = NULL, smoothed=FALSE)
                        {
# reshape
df <- as.data.frame(twisst_object$weights[[chr]])

if (smoothed==TRUE){
  names(twisst_object$pos)<-names(twisst_object$weights)
}
  
df$pos <- unlist(twisst_object$pos[[chr]])

df_long <- df %>%
  pivot_longer(-pos, names_to = "topo", values_to = "weight")

# basic overlaid areas
ggplot(df_long, aes(x = pos/1e6, y = weight, fill = topo)) +
  geom_area(alpha = 0.6, position = "identity") +
  scale_x_continuous(expand = c(0.01, 0), name = "Chromosome Position (Mb)",  breaks = scales::breaks_extended(n = 4), limits=xlim/1e6)+
  labs(x = "Position", y = "Topology weight") +
  ylim(0, 1)+
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
}

plot.twisst.gg.topo.subset <- function(twisst_object,chr,
                           xlim = NULL,topo, smoothed=FALSE)
{
  # reshape
  df <- as.data.frame(twisst_object$weights[[chr]][[topo]])
  
  if (smoothed==TRUE){
    names(twisst_object$pos)<-names(twisst_object$weights)
  }
  
  df$pos <- unlist(twisst_object$pos[[chr]])
  
  df_long <- df %>%
    pivot_longer(-pos, names_to = "topo", values_to = "weight")
  
  # basic overlaid areas
  ggplot(df_long, aes(x = pos/1e6, y = weight, fill = topo)) +
    geom_area(alpha = 0.6, position = "identity") +
    scale_x_continuous(expand = c(0.01, 0), name = "Chromosome Position (Mb)",  breaks = scales::breaks_extended(n = 4), limits=xlim/1e6)+
    labs(x = "Position", y = "Topology weight") +
    ylim(0, 1)+
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
      legend.position = "none",
      legend.text = element_text(size=rel(1.2), family = "Arial"),
      panel.spacing = unit(2, "lines")
    ) +
    labs(y = expression("Topology weight"))
}