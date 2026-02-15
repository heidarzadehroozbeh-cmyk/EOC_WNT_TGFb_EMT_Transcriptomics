 # 05_Volcano Plots.R
library(ggplot2)
library(ggrepel)

# Parameters
GSE_list <- c("GSE14407", "GSE38666", "GSE52037", "GSE216150")  # Replace #### with your 4th dataset
FDR_thresh <- 0.05
logFC_thresh <- 1

# Ensure results/plots exists
plots_dir <- file.path("results", "plots")
if(!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# Loop over all datasets
for (gse_id in GSE_list) {
  
  deg_file <- file.path("results/tables", paste0(gse_id, "_limma_all.tsv"))
  plot_file <- file.path(plots_dir, paste0(gse_id, "_volcano_labeled.png"))
  
  # Load limma results
  deg <- read.table(deg_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Volcano plot colors
  deg$color <- "grey"
  deg$color[deg$adj.P.Val < FDR_thresh & deg$logFC > logFC_thresh] <- "red"
  deg$color[deg$adj.P.Val < FDR_thresh & deg$logFC < -logFC_thresh] <- "blue"
  
  # Label top 10 genes by smallest FDR
  top_genes <- deg[deg$adj.P.Val < FDR_thresh, ]
  top_genes <- top_genes[order(top_genes$adj.P.Val), ]
  top_genes <- head(top_genes, 10)
  
  # Plot
  volcano_plot <- ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = color), alpha = 0.6) +
    scale_color_identity() +
    geom_text_repel(data = top_genes, aes(label = gene), max.overlaps = 20) +
    theme_minimal() +
    labs(
      title = paste("Volcano plot:", gse_id),
      x = "log2 Fold Change",
      y = "-log10 Adjusted P-value"
    )
  
  # Save plot
  ggsave(plot_file, volcano_plot, width = 7, height = 5)
  message("Volcano plot saved to: ", plot_file)
}
