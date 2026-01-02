#!/usr/bin/env Rscript

#' Standalone Test: DESandbox with Airway Dataset
#' 
#' This script tests DESandbox functionality by sourcing the R files directly
#' without requiring package installation. Perfect for development testing.

cat("=== DESandbox Standalone Test with Airway Dataset ===\n\n")

# ========================================
# Load Required Libraries
# ========================================
cat("Loading required libraries...\n")

required_packages <- c(
  "DESeq2", "edgeR", "limma", 
  "SummarizedExperiment", "S4Vectors",
  "airway", "methods"
)

missing_packages <- c()
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing_packages <- c(missing_packages, pkg)
  }
}

if (length(missing_packages) > 0) {
  cat("\nMissing required packages:\n")
  cat(paste("  -", missing_packages, collapse = "\n"), "\n\n")
  cat("Install with:\n")
  cat("  BiocManager::install(c(", paste0("'", missing_packages, "'", collapse = ", "), "))\n\n")
  stop("Please install missing packages first.")
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(airway)
  library(methods)
})

cat("✓ All required packages loaded\n\n")

# ========================================
# Source DESandbox Functions
# ========================================
cat("Loading DESandbox functions from R/ directory...\n")

r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
if (length(r_files) == 0) {
  stop("Cannot find R/ directory. Please run from package root directory.")
}

for (f in r_files) {
  source(f)
  cat(sprintf("  ✓ Loaded %s\n", basename(f)))
}

cat("\n")

# ========================================
# Load and Prepare Data
# ========================================
cat("Step 1: Loading airway dataset...\n")
data(airway)

counts <- assay(airway)
metadata <- as.data.frame(colData(airway))
metadata$condition <- metadata$dex
rownames(metadata) <- metadata$Run

# Filter for faster test
counts <- counts[rowSums(counts) > 100, ]

cat(sprintf("  - Loaded %d genes across %d samples\n", nrow(counts), ncol(counts)))
cat(sprintf("  - Conditions: %s\n", paste(unique(metadata$condition), collapse = ", ")))

# ========================================
# Create DESandbox Object
# ========================================
cat("\nStep 2: Creating DESandbox object...\n")
dso <- create_desandbox_object(
  counts = counts,
  metadata = metadata,
  condition_column = "condition"
)

cat(sprintf("  ✓ Created SummarizedExperiment with %d genes and %d samples\n", 
            nrow(dso), ncol(dso)))

# ========================================
# Run DESeq2
# ========================================
cat("\nStep 3: Running DESeq2 analysis...\n")
deseq2_results <- run_deseq2(dso, alpha = 0.05)

cat(sprintf("  ✓ DESeq2: Analyzed %d genes\n", nrow(deseq2_results$results)))

# ========================================
# Run edgeR
# ========================================
cat("\nStep 4: Running edgeR analysis...\n")
edger_results <- run_edger(dso, normalization = "TMM")

cat(sprintf("  ✓ edgeR: Analyzed %d genes\n", nrow(edger_results$results)))

# ========================================
# Run limma-voom
# ========================================
cat("\nStep 5: Running limma-voom analysis...\n")
limma_results <- run_limma_voom(dso, normalization = "TMM")

cat(sprintf("  ✓ limma-voom: Analyzed %d genes\n", nrow(limma_results$results)))

# ========================================
# Combine Results
# ========================================
results <- list(
  DESeq2 = deseq2_results,
  edgeR = edger_results,
  `limma-voom` = limma_results
)

# ========================================
# Filter Significant Genes
# ========================================
cat("\nStep 6: Filtering significant genes (padj < 0.05, |log2FC| > 1)...\n")

filtered_results <- lapply(results, function(res) {
  filter_results(
    res$results,
    padj_threshold = 0.05,
    lfc_threshold = 1
  )
})

for (method in names(filtered_results)) {
  n_sig <- nrow(filtered_results[[method]])
  cat(sprintf("  - %s: %d significant genes\n", method, n_sig))
}

# ========================================
# Compare Methods
# ========================================
cat("\nStep 7: Comparing methods...\n")

comparison <- compare_methods(
  results,
  padj_threshold = 0.05,
  lfc_threshold = 1
)

cat("\n  Summary Statistics:\n")
print(comparison$summary_stats)

cat(sprintf("\n  Consensus genes (all methods): %d\n", 
            length(comparison$consensus_genes)))

# ========================================
# Create Consensus List
# ========================================
cat("\nStep 8: Creating consensus gene list (min 2 methods)...\n")

consensus <- create_consensus_list(
  results,
  min_methods = 2,
  padj_threshold = 0.05,
  lfc_threshold = 1
)

cat(sprintf("  ✓ Found %d consensus genes\n", nrow(consensus)))

if (nrow(consensus) > 0) {
  cat("\n  Top 5 Consensus Genes:\n")
  print(head(consensus[, c("gene", "n_methods_support", "log2FoldChange", "min_padj")], 5))
}

# ========================================
# Export Results
# ========================================
cat("\nStep 9: Exporting results...\n")

output_dir <- "test_output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

export_results(
  results,
  output_prefix = file.path(output_dir, "test_analysis"),
  formats = c("csv", "json", "rds")
)

write.csv(
  consensus,
  file = file.path(output_dir, "consensus_genes.csv"),
  row.names = FALSE
)

cat(sprintf("  ✓ Results exported to '%s/'\n", output_dir))

# ========================================
# Test Caching
# ========================================
cat("\nStep 10: Testing caching functionality...\n")

cache_dir <- file.path(output_dir, "cache")
cache_path <- cache_analysis(
  results,
  cache_dir = cache_dir,
  cache_name = "test_cache"
)

cat(sprintf("  ✓ Cached to: %s\n", basename(cache_path)))

# Load cached results
cached <- load_cached_analysis(cache_path, check_version = FALSE)
cat("  ✓ Successfully loaded cached results\n")

# ========================================
# Summary
# ========================================
cat("\n=== Test Complete! ===\n\n")
cat("Summary:\n")
cat(sprintf("  ✓ All 3 DE methods ran successfully\n"))
cat(sprintf("  ✓ Found %d consensus genes\n", nrow(consensus)))
cat(sprintf("  ✓ Results exported to %s/\n", output_dir))
cat(sprintf("  ✓ Caching system working\n"))

cat("\nOutput Files:\n")
output_files <- list.files(output_dir, recursive = TRUE)
for (f in output_files) {
  cat(sprintf("  - %s\n", f))
}

cat("\n=== All Tests Passed! ===\n")
