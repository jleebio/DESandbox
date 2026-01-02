#!/usr/bin/env Rscript

# DESandbox Visualization Examples
# Demonstrates volcano plots, Venn diagrams, and comparison plots

# Load required packages
library(airway)
library(SummarizedExperiment)

# Source DESandbox functions
source("R/data_validation.R")
source("R/de_analysis.R")
source("R/result_processing.R")
source("R/comparison.R")
source("R/visualization.R")

cat("=== DESandbox Visualization Demo ===\n\n")

# ============================================================================
# 1. Load Data and Run Analysis
# ============================================================================
cat("Step 1: Loading airway dataset and running DE analysis...\n")

data(airway)
counts_matrix <- assay(airway, "counts")
metadata_df <- as.data.frame(colData(airway))

dso <- create_desandbox_object(
  counts = counts_matrix,
  metadata = metadata_df,
  condition_column = "dex"
)

# Run all three methods
cat("Running DESeq2, edgeR, and limma-voom...\n")
all_results <- run_all_methods(dso)

cat(sprintf("✓ Analysis complete\n\n"))

# ============================================================================
# 2. Prepare Data for Visualization
# ============================================================================
cat("Step 2: Preparing unified results data frame...\n")

# Combine results from all methods into a single data frame
combined_df <- data.frame(rbind(
  data.frame(
    method = "DESeq2",
    gene_id = all_results$DESeq2$results$gene,
    log2FC = all_results$DESeq2$results$log2FoldChange,
    padj = all_results$DESeq2$results$padj
  ),
  data.frame(
    method = "edgeR",
    gene_id = all_results$edgeR$results$gene,
    log2FC = all_results$edgeR$results$log2FoldChange,
    padj = all_results$edgeR$results$padj
  ),
  data.frame(
    method = "limma-voom",
    gene_id = all_results$`limma-voom`$results$gene,
    log2FC = all_results$`limma-voom`$results$log2FoldChange,
    padj = all_results$`limma-voom`$results$padj
  )
))

cat(sprintf("✓ Combined data frame: %d rows, %d methods\n\n", 
            nrow(combined_df), 
            length(unique(combined_df$method))))

# ============================================================================
# 3. Create Volcano Plots
# ============================================================================
cat("Step 3: Creating volcano plots...\n")

# Create output directory
if (!dir.exists("visualization_output")) {
  dir.create("visualization_output")
}

# 3a. Faceted volcano plot (all methods)
cat("  Creating faceted volcano plot...\n")
p_volcano_facet <- plot_volcano(
  combined_df,
  log2FC_threshold = 1,
  padj_threshold = 0.05,
  facet_by_method = TRUE,
  label_top_n = 5,
  title = "Differential Expression: Airway Dataset",
  subtitle = "Dexamethasone treatment vs. control"
)

ggplot2::ggsave(
  "visualization_output/volcano_faceted.png",
  p_volcano_facet,
  width = 12,
  height = 4,
  dpi = 300
)
cat("  ✓ Saved: visualization_output/volcano_faceted.png\n")

# 3b. Single volcano plot (DESeq2 only)
cat("  Creating DESeq2 volcano plot...\n")
deseq2_df <- combined_df[combined_df$method == "DESeq2", ]

p_volcano_deseq2 <- plot_volcano(
  deseq2_df,
  log2FC_threshold = 1.5,
  padj_threshold = 0.01,
  facet_by_method = FALSE,
  label_top_n = 15,
  label_by = "padj",
  title = "DESeq2 Analysis",
  subtitle = "Stringent thresholds: |log2FC| > 1.5, padj < 0.01"
)

ggplot2::ggsave(
  "visualization_output/volcano_deseq2.png",
  p_volcano_deseq2,
  width = 8,
  height = 6,
  dpi = 300
)
cat("  ✓ Saved: visualization_output/volcano_deseq2.png\n\n")

# ============================================================================
# 4. Create Venn Diagrams
# ============================================================================
cat("Step 4: Creating Venn diagrams...\n")

# 4a. Venn diagram for all DE genes (both directions)
cat("  Creating Venn diagram (all DE genes)...\n")
p_venn_all <- plot_venn_de(
  combined_df,
  log2FC_threshold = 1,
  padj_threshold = 0.05,
  direction = "both",
  title = "DE Gene Overlap: All Significant Genes"
)

ggplot2::ggsave(
  "visualization_output/venn_all_de.png",
  p_venn_all,
  width = 8,
  height = 8,
  dpi = 300
)
cat("  ✓ Saved: visualization_output/venn_all_de.png\n")

# 4b. Venn diagram for upregulated genes only
cat("  Creating Venn diagram (upregulated genes)...\n")
p_venn_up <- plot_venn_de(
  combined_df,
  log2FC_threshold = 1,
  padj_threshold = 0.05,
  direction = "up",
  title = "DE Gene Overlap: Upregulated Genes"
)

ggplot2::ggsave(
  "visualization_output/venn_upregulated.png",
  p_venn_up,
  width = 8,
  height = 8,
  dpi = 300
)
cat("  ✓ Saved: visualization_output/venn_upregulated.png\n")

# 4c. Venn diagram for downregulated genes only
cat("  Creating Venn diagram (downregulated genes)...\n")
p_venn_down <- plot_venn_de(
  combined_df,
  log2FC_threshold = 1,
  padj_threshold = 0.05,
  direction = "down",
  title = "DE Gene Overlap: Downregulated Genes"
)

ggplot2::ggsave(
  "visualization_output/venn_downregulated.png",
  p_venn_down,
  width = 8,
  height = 8,
  dpi = 300
)
cat("  ✓ Saved: visualization_output/venn_downregulated.png\n\n")

# ============================================================================
# 5. Get Gene Sets from Venn Analysis
# ============================================================================
cat("Step 5: Extracting gene sets from Venn analysis...\n")

# Get actual gene lists
gene_sets_up <- plot_venn_de(
  combined_df,
  direction = "up",
  return_type = "sets"
)

gene_sets_down <- plot_venn_de(
  combined_df,
  direction = "down",
  return_type = "sets"
)

cat("Upregulated gene counts by method:\n")
for (method in names(gene_sets_up)) {
  cat(sprintf("  %s: %d genes\n", method, length(gene_sets_up[[method]])))
}

cat("\nDownregulated gene counts by method:\n")
for (method in names(gene_sets_down)) {
  cat(sprintf("  %s: %d genes\n", method, length(gene_sets_down[[method]])))
}

# Find consensus genes (all 3 methods agree)
consensus_up <- Reduce(intersect, gene_sets_up)
consensus_down <- Reduce(intersect, gene_sets_down)

cat(sprintf("\n✓ Consensus upregulated genes (all 3 methods): %d\n", 
            length(consensus_up)))
cat(sprintf("✓ Consensus downregulated genes (all 3 methods): %d\n\n", 
            length(consensus_down)))

# ============================================================================
# 6. Create DE Gene Count Comparison Plot
# ============================================================================
cat("Step 6: Creating DE gene count comparison plot...\n")

p_counts <- plot_de_counts(
  combined_df,
  log2FC_threshold = 1,
  padj_threshold = 0.05,
  split_by_direction = TRUE
)

ggplot2::ggsave(
  "visualization_output/de_counts_comparison.png",
  p_counts,
  width = 8,
  height = 6,
  dpi = 300
)
cat("  ✓ Saved: visualization_output/de_counts_comparison.png\n\n")

# ============================================================================
# 7. Advanced: Custom Volcano Plot Styling
# ============================================================================
cat("Step 7: Creating custom-styled volcano plot...\n")

# Custom colors
custom_colors <- c(
  "Upregulated" = "#D62728",      # Darker red
  "Downregulated" = "#1F77B4",    # Darker blue
  "Not significant" = "#D3D3D3"   # Light gray
)

p_volcano_custom <- plot_volcano(
  deseq2_df,
  log2FC_threshold = 1,
  padj_threshold = 0.05,
  facet_by_method = FALSE,
  label_top_n = 20,
  point_size = 2,
  point_alpha = 0.7,
  colors = custom_colors,
  title = "Custom Styled Volcano Plot",
  subtitle = "DESeq2 results with custom aesthetics"
)

ggplot2::ggsave(
  "visualization_output/volcano_custom.png",
  p_volcano_custom,
  width = 10,
  height = 7,
  dpi = 300
)
cat("  ✓ Saved: visualization_output/volcano_custom.png\n\n")

# ============================================================================
# 8. Export Gene Lists
# ============================================================================
cat("Step 8: Exporting gene lists...\n")

# Export consensus genes
write.csv(
  data.frame(gene_id = consensus_up),
  "visualization_output/consensus_upregulated.csv",
  row.names = FALSE
)

write.csv(
  data.frame(gene_id = consensus_down),
  "visualization_output/consensus_downregulated.csv",
  row.names = FALSE
)

cat("  ✓ Saved: visualization_output/consensus_upregulated.csv\n")
cat("  ✓ Saved: visualization_output/consensus_downregulated.csv\n\n")

# ============================================================================
# Summary
# ============================================================================
cat("=== Visualization Demo Complete ===\n\n")
cat("Generated visualizations:\n")
cat("  1. volcano_faceted.png      - Faceted volcano plot (all methods)\n")
cat("  2. volcano_deseq2.png        - DESeq2 volcano plot with labels\n")
cat("  3. venn_all_de.png           - Venn diagram (all DE genes)\n")
cat("  4. venn_upregulated.png      - Venn diagram (upregulated)\n")
cat("  5. venn_downregulated.png    - Venn diagram (downregulated)\n")
cat("  6. de_counts_comparison.png  - Bar plot of DE gene counts\n")
cat("  7. volcano_custom.png        - Custom styled volcano plot\n\n")

cat("Exported gene lists:\n")
cat("  - consensus_upregulated.csv\n")
cat("  - consensus_downregulated.csv\n\n")

cat("All files saved to: visualization_output/\n")
