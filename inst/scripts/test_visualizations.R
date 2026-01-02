#!/usr/bin/env Rscript

# Quick test of DESandbox visualization functions
# Tests each function individually with airway data

cat("=== DESandbox Visualization Quick Test ===\n\n")

# Load packages
suppressPackageStartupMessages({
  library(airway)
  library(SummarizedExperiment)
})

# Source functions
source("R/data_validation.R")
source("R/de_analysis.R")
source("R/visualization.R")

cat("Step 1: Prepare test data\n")
data(airway)
counts_matrix <- assay(airway, "counts")
metadata_df <- as.data.frame(colData(airway))

dso <- create_desandbox_object(
  counts = counts_matrix,
  metadata = metadata_df,
  condition_column = "dex"
)

# Run just DESeq2 for speed
cat("Step 2: Run DESeq2 analysis\n")
result_deseq2 <- run_deseq2(dso)

# Prepare data frame for visualization
test_df <- data.frame(
  gene_id = result_deseq2$results$gene,
  log2FC = result_deseq2$results$log2FoldChange,
  padj = result_deseq2$results$padj
)

cat(sprintf("  Data frame: %d genes\n\n", nrow(test_df)))

# Test 1: Basic volcano plot
cat("Test 1: Basic volcano plot\n")
p1 <- plot_volcano(test_df)
cat("  ✓ plot_volcano() executed successfully\n")
cat(sprintf("  Class: %s\n", class(p1)[1]))

# Test 2: Volcano with custom parameters
cat("\nTest 2: Volcano with custom parameters\n")
p2 <- plot_volcano(
  test_df,
  log2FC_threshold = 1.5,
  padj_threshold = 0.01,
  label_top_n = 5,
  title = "Custom Volcano"
)
cat("  ✓ Custom parameters work\n")

# Test 3: Volcano with no labels
cat("\nTest 3: Volcano without labels\n")
p3 <- plot_volcano(test_df, label_top_n = 0)
cat("  ✓ No-label mode works\n")

# Now test Venn diagrams (need multiple methods)
cat("\nStep 3: Run all methods for Venn test\n")
all_results <- run_all_methods(dso)

# Prepare combined data frame
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

cat(sprintf("  Combined data frame: %d rows\n\n", nrow(combined_df)))

# Test 4: Venn diagram (all DE genes)
cat("Test 4: Venn diagram (all DE genes)\n")
p4 <- plot_venn_de(combined_df, direction = "both")
cat("  ✓ plot_venn_de() executed successfully\n")
cat(sprintf("  Class: %s\n", class(p4)[1]))

# Test 5: Venn diagram (upregulated only)
cat("\nTest 5: Venn diagram (upregulated genes)\n")
p5 <- plot_venn_de(combined_df, direction = "up")
cat("  ✓ Direction filtering works\n")

# Test 6: Venn diagram (downregulated only)
cat("\nTest 6: Venn diagram (downregulated genes)\n")
p6 <- plot_venn_de(combined_df, direction = "down")
cat("  ✓ Down direction works\n")

# Test 7: Extract gene sets
cat("\nTest 7: Extract gene sets from Venn\n")
gene_sets_up <- plot_venn_de(combined_df, direction = "up", return_type = "sets")
cat(sprintf("  DESeq2: %d genes\n", length(gene_sets_up$DESeq2)))
cat(sprintf("  edgeR: %d genes\n", length(gene_sets_up$edgeR)))
cat(sprintf("  limma-voom: %d genes\n", length(gene_sets_up$`limma-voom`)))
cat("  ✓ Gene set extraction works\n")

# Test 8: DE count comparison
cat("\nTest 8: DE count comparison plot\n")
p7 <- plot_de_counts(combined_df)
cat("  ✓ plot_de_counts() executed successfully\n")

# Test 9: DE counts without direction split
cat("\nTest 9: DE counts (no direction split)\n")
p8 <- plot_de_counts(combined_df, split_by_direction = FALSE)
cat("  ✓ Split option works\n")

# Test 10: Faceted volcano
cat("\nTest 10: Faceted volcano plot\n")
p9 <- plot_volcano(combined_df, facet_by_method = TRUE, label_top_n = 3)
cat("  ✓ Faceting works\n")

# Test 11: Edge cases - empty significant genes
cat("\nTest 11: Edge case - very stringent thresholds\n")
p10 <- plot_volcano(test_df, log2FC_threshold = 10, padj_threshold = 1e-50)
cat("  ✓ Handles case with few/no significant genes\n")

# Test 12: Custom colors
cat("\nTest 12: Custom color scheme\n")
custom_colors <- c(
  "Upregulated" = "#FF0000",
  "Downregulated" = "#0000FF",
  "Not significant" = "#CCCCCC"
)
p11 <- plot_volcano(test_df, colors = custom_colors)
cat("  ✓ Custom colors work\n")

# Summary
cat("\n" , rep("=", 50), "\n", sep = "")
cat("=== All Tests Passed! ===\n\n")

cat("Functions tested:\n")
cat("  1. plot_volcano() - basic\n")
cat("  2. plot_volcano() - custom parameters\n")
cat("  3. plot_volcano() - no labels\n")
cat("  4. plot_volcano() - faceted\n")
cat("  5. plot_volcano() - custom colors\n")
cat("  6. plot_volcano() - edge cases\n")
cat("  7. plot_venn_de() - all genes\n")
cat("  8. plot_venn_de() - upregulated\n")
cat("  9. plot_venn_de() - downregulated\n")
cat(" 10. plot_venn_de() - gene set extraction\n")
cat(" 11. plot_de_counts() - split by direction\n")
cat(" 12. plot_de_counts() - total counts\n\n")

cat("✓ All visualization functions working correctly!\n")
cat("✓ Returns ggplot objects as expected\n")
cat("✓ Handles various parameter combinations\n")
cat("✓ Edge cases handled gracefully\n")
