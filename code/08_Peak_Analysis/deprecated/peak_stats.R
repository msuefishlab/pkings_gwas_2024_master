library(dplyr)
library(FSA) # For Dunn's test

# Function to analyze each peak separately
analyze_peak <- function(peak_name, data) {
  cat("\n### Analyzing", peak_name, "###\n")
  
  df_peak <- data %>% filter(peak == peak_name, category == "peak")
  
  # Step 1: Compare BP vs TP across all populations in peak regions
  test1 <- wilcox.test(avg_pi ~ phenotype_pop1, data = df_peak, alternative = "less")  # Test for BP < TP
  cat("Step 1 - Wilcoxon Test (BP < TP across populations in peak regions): p =", test1$p.value, "\n")
  
  significant_bp_populations <- c()  # Track BP populations with significantly reduced π
  
  if (test1$p.value < 0.05) {
    cat("Significant reduction in BP vs TP found. Checking which BP populations contribute...\n")
    
    # Step 2: Test for variation across populations
    test2 <- kruskal.test(avg_pi ~ population, data = df_peak)
    cat("Step 2 - Kruskal-Wallis Test (Variation among populations in peak regions): p =", test2$p.value, "\n")
    
    if (test2$p.value < 0.05) {
      cat("Performing post-hoc Dunn's test to identify significant BP populations...\n")
      dunn_results <- dunnTest(avg_pi ~ population, data = df_peak, method="bonferroni")
      print(dunn_results)
      
      # Extract population pairs involved in significant comparisons
      significant_pairs <- dunn_results$res %>%
        filter(P.adj < 0.05) %>%
        pull(Comparison)
      
      # Extract population names from comparisons
      population_pairs <- strsplit(significant_pairs, " - ")
      
      # Determine which BP populations have significantly lower π than TP populations
      for (pair in population_pairs) {
        pop1 <- pair[1]
        pop2 <- pair[2]
        
        # Get mean π for each population
        pi_means <- df_peak %>%
          filter(population %in% c(pop1, pop2)) %>%
          group_by(population) %>%
          summarize(mean_pi = mean(avg_pi), phenotype = unique(phenotype_pop1), .groups = "drop")
        
        # Ensure BP has lower π than TP
        bp_pop <- pi_means %>% filter(phenotype == "BP")
        tp_pop <- pi_means %>% filter(phenotype == "TP")
        
        if (nrow(bp_pop) == 1 && nrow(tp_pop) == 1 && bp_pop$mean_pi < tp_pop$mean_pi) {
          significant_bp_populations <- c(significant_bp_populations, bp_pop$population)
        }
      }
      
      significant_bp_populations <- unique(significant_bp_populations)
      cat("BP populations where π is significantly reduced:", paste(significant_bp_populations, collapse=", "), "\n")
    }
  }
  
  # Step 3: Compare peak vs. random within significant BP populations
  if (length(significant_bp_populations) > 0) {
    cat("\nStep 3 - Checking peak vs. random regions for significant BP populations...\n")
    
    for (pop in significant_bp_populations) {
      df_pop <- data %>% filter(peak == peak_name, population == pop)
      
      if (length(unique(df_pop$category)) == 2) {  # Ensure both peak and random exist
        test3 <- wilcox.test(avg_pi ~ category, data = df_pop, alternative = "less")  # Test for Peak < Random
        cat("\nStep 3 - Wilcoxon Test (Peak < Random in", pop, "): p =", test3$p.value, "\n")
      }
    }
  } else {
    cat("\nNo BP populations showing significant π reduction. Skipping further analysis.\n")
  }
}
