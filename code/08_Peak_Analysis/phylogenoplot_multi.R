library(ggtree)
library(phangorn)
library(tidyverse)

create_phylo_geno_plot_multi<- function(tree,snpdata,peaks_per_plot = 4)
{
  old<-tree$tip.label
  new<-sapply(old,function(x) strsplit(x, "_")[[1]][1])
  names(new)<-NULL
  
  nameframe<-data.frame(old.labels=old,new.labels=new)
  
  tree$tip.label[tree$tip.label %in% nameframe$old.labels] <- nameframe$new.labels
  
  
  ### recover pop tree
  
  name_to_pop<-snpdata %>% dplyr::select(name,population) %>% unique()
  name_to_pop <- name_to_pop %>% mutate(population = as.character(population))
  
  name_to_pop <- name_to_pop %>%
    mutate(tip_index = match(population, tree$tip.label))
  
  # Function to create a polytomy for individuals in a given population
  add_polytomy_to_population <- function(tree, individuals, population_label) {
    # Create a star tree (polytomy) with individuals
    polytomy <- stree(length(individuals), tip.label = individuals)
    
    # Find the current index of the population tip by its label
    parent_tip_index <- which(tree$tip.label == population_label)
    
    # Add the polytomy to the tree
    new_tree <- bind.tree(tree, polytomy, where = parent_tip_index, position = 0)
    
    return(new_tree)
  }
  
  # Add polytomies for each population
  for (pop in unique(name_to_pop$population)) {
    # Get individuals for this population
    individuals <- name_to_pop %>%
      filter(population == pop) %>%
      pull(name)
    
    # Add the polytomy to the tree
    tree <- add_polytomy_to_population(tree, individuals, pop)
  }
  
  test_data<-as.data.frame(snpdata)
  pheno_data<-test_data[c(1,9,10)]
  
  
  
  # Create a data frame to identify tips that should be hidden
  tip_data <- data.frame(
    label = tree$tip.label,
    is_polytomy = tree$tip.label %in% name_to_pop$name # Identify tips belonging to polytomies
  )
  
  # Identify ancestor nodes for each population
  ancestor_nodes <- name_to_pop %>%
    group_by(population) %>%
    summarize(ancestor_node = getMRCA(tree, name)) %>%
    mutate(label = population)
  
  ordered_levels <- test_data$peak[order(test_data$peak)]
  test_data$peak <- factor(test_data$peak, levels = unique(ordered_levels))
  
  peak_levels <- levels(test_data$peak)
  peak_groups <- split(peak_levels, ceiling(seq_along(peak_levels) / peaks_per_plot))
  plot_list <- list()
  
  for (group in peak_groups) {
  # Plot the tree
  p <- ggtree(tree, color = "white", size = 1)
  
  
  p <- p %<+% pheno_data + geom_tippoint(aes(fill=Phenotype),shape = 22, size = 1, color = "black",stroke=0.02,  position = position_nudge(x = 3))
  
  # Extract population boundaries for all peaks
  population_boundaries <- p$data %>%
    filter(!is.na(population)) %>%
    #group_by(peak, population) %>%
    summarize(y_min = min(y), y_max = max(y), .groups = 'drop') %>%
    mutate(boundary = y_max + 0.5)
  
  # Create heatmaps for each peak
  for (i in group) {
      filtered_test_data<- test_data %>% dplyr::select(-Phenotype,-population,-maxp) %>% filter(peak==i)
      nested_test_data<-NULL
      nested_test_data<-nest(filtered_test_data,pos=pos,chr=chr,geno=geno,snp_order_within_individual=snp_order_within_individual)
      p <- p %<+% nested_test_data
      p<-  p+ geom_facet(
          panel = paste0("Peak ", i),
          data = td_unnest(c(snp_order_within_individual, geno)),
          geom = geom_tile,
          mapping = aes(x = snp_order_within_individual, fill = geno)
        ) +
        scale_fill_manual(
          name = "Legend",
          values = c("NA" = "black", "0" = "#cccccc", "1" = "#ffce95", "2" = "#fd7660", "BP" = "blue", "TP" = "red", "Wobble" = "purple"),
          labels = c("NA" = "Missing", "0" = "REF", "1" = "HET", "2" = "ALT", "BP" = "Biphasic", "TP" = "Triphasic", "Wobble" = "Intermediate")
        ) 
      print(i)
      p$data<-p$data %>% dplyr::select(-peak,-pos,-chr,-geno,-snp_order_within_individual)
  }
  
  p<-p+
    geom_hline(
      data = population_boundaries %>% mutate(.panel = paste0("Peak ", i)),
      aes(yintercept = boundary),
      color = "black"
    )
  # Extract ggtree plot data and merge with tip_data
  p$data <- p$data %>%
    left_join(tip_data, by = "label") %>%
    mutate(is_polytomy = ifelse(is.na(is_polytomy), FALSE, is_polytomy))
  
  # Filter out polytomy tips and branches for plotting
  non_polytomy_data <- p$data %>% filter(!is_polytomy)
  
  # Add population labels to ancestor nodes
  p <- p + geom_text2(
    data = p$data %>% filter(node %in% ancestor_nodes$ancestor_node),
    aes(label = ancestor_nodes$label[match(node, ancestor_nodes$ancestor_node)]),
    hjust = -0.3,
    size = 4
  )
  
  # Plot visible branches (non-polytomy)
  p <- p + geom_tree(data = non_polytomy_data)+theme_tree2()+
    coord_cartesian(clip = "off") + 
    theme(strip.background = element_blank(), 
          legend.position = "right",          # Move legend to the bottom
          legend.direction = "horizontal",     # Arrange legend items horizontally
          legend.box = "horizontal",          # Ensure the legend box aligns horizontally)
          plot.margin = margin(10, 10, 10, 10)
    )
  plot_list <- c(plot_list, list(p))
  }
  

  return(plot_list)
}

