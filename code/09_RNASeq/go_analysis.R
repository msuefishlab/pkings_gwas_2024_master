library(topGO)


run_go_analysis <- function(de_genes, gene_to_go, ontology = "CC", p_threshold = 0.05) {
  # Create a named binary vector: 1 if gene is in de_genes, 0 otherwise.
  gene_list <- factor(as.integer(names(gene_to_go) %in% de_genes))
  names(gene_list) <- names(gene_to_go)
  
  # Create the topGOdata object.
  go_data <- new(
    "topGOdata",
    ontology = ontology,
    allGenes = gene_list,
    geneSelectionFun = function(x) x == 1,
    annot = annFUN.gene2GO,
    gene2GO = gene_to_go
  )
  
  # Run the GO enrichment test.
  go_results <- runTest(go_data, algorithm = "classic", statistic = "fisher")
  
  # Generate the results table for all nodes.
  all_results <- GenTable(
    go_data,
    classicFisher = go_results,
    orderBy = "classicFisher",
    ranksOf = "classicFisher",
    topNodes = length(score(go_results))
  )
  
  # Convert p-values to numeric.
  all_results$classicFisher <- as.numeric(all_results$classicFisher)
  
  # Filter for significant GO terms.
  filtered_results <- all_results[all_results$classicFisher < p_threshold, ]
  
  return(filtered_results)
}

plot_go_results <- function(filtered_results, plot_title = "GO Enrichment Bar Plot") {
  ggplot(filtered_results, aes(x = reorder(Term, -log10(classicFisher)), y = -log10(classicFisher))) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = Significant), hjust = -0.2, size = 3) +
    coord_flip() +
    labs(
      title = plot_title,
      x = "GO Terms",
      y = "-log10(P-value)"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 14)
    )
}

plot_combined_go_results <- function(results_up, results_down ) {
  # Add a direction label to each result set
  results_up$direction <- "TP"
  results_down$direction <- "BP"
  
  # Create a signed enrichment score:
  # For upregulated genes, enrichment = -log10(p)
  results_up$signed_logP <- -log10(results_up$classicFisher)
  # For downregulated genes, enrichment = -(-log10(p)) = log10(p) * (-1)
  results_down$signed_logP <- -(-log10(results_down$classicFisher))
  # (Alternatively, you can write:
  # results_down$signed_logP <- -log10(results_down$classicFisher) * -1
  # which yields the same result.)
  
  # Combine the two data frames
  combined_results <- rbind(results_up, results_down)
  
  # Optionally, reorder the GO terms by the signed enrichment value
  # so that the most extreme (either positive or negative) are at the top.
  # Here, we order by the absolute value of signed_logP in decreasing order.
  ordered_levels <- unique(
    combined_results[order(combined_results$signed_logP, decreasing = TRUE), "Term"]
  )
  combined_results$Term <- factor(combined_results$Term, levels = ordered_levels)
  
  # Create the plot
  ggplot(combined_results, aes(x = Term, y = signed_logP, fill = direction)) +
    geom_col() +
    coord_flip() +  # horizontal bars
    labs(
      x = "GO Terms",
      y = "-log10(p-value)"
    ) +
    geom_text(aes(
      label = Significant,
      hjust = ifelse(signed_logP > 0, -0.2, 1.2)  # adjust hjust based on the sign
    ), size = 3) +    theme_minimal() +
    theme_minimal() +
    theme(
      legend.position = "none",  # Hide legend if needed
      panel.grid.major.y = element_blank(),  # Remove major vertical grid lines
      panel.grid.minor.y = element_blank(),  # Remove major vertical grid lines
      panel.grid.major.x = element_blank(),  # Remove major vertical grid lines
      panel.grid.minor.x = element_blank(),  # Remove minor vertical grid lines
      axis.title.y=element_blank(),
      axis.line.x = element_line(color="black", size = 0.5),
      axis.text.y = element_text(size = 10),
    ) +
    scale_fill_manual(values = c("TP" = "red", "BP" = "blue"))
}

