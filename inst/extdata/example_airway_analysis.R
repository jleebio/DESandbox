#!/usr/bin/env Rscript

#' Complete Working Example: DESandbox with Airway Dataset
#' 
#' This script demonstrates a complete differential expression analysis workflow
#' using the airway dataset from Bioconductor.
#' 
#' Dataset: RNA-seq experiment studying the effect of dexamethasone treatment
#' on airway smooth muscle cells.

# Load required libraries
library(DESandbox)
library(airway)

cat("=== DESandbox Example Analysis ===\n\n")

# ========================================
# Step 1: Load and Prepare Data
# ========================================
cat("Step 1: Loading airway dataset...\n")
data(airway)

# Extract counts and metadata
counts <- assay(airway)
metadata <- as.data.frame(colData(airway))

# Prepare metadata for analysis
metadata$condition <- metadata$dex
rownames(metadata) <- metadata$Run

# Filter low-count genes
counts <- counts[rowSums(counts) > 100, ]

cat(sprintf("  - Loaded %d genes across %d samples\n", nrow(counts), ncol(counts)))
cat(sprintf("  - Conditions: %s\n", paste(levels(factor(metadata$condition)), collapse = ", ")))

# ========================================
# Step 2: Create DESandbox Object
# ========================================
cat("\nStep 2: Creating DESandbox object...\n")
dso <- create_desandbox_object(
  counts = counts,
  metadata = metadata,
  condition_column = "condition",
  validate = TRUE
)

cat(sprintf("  - Created SummarizedExperiment with %d genes and %d samples\n", 
            nrow(dso), ncol(dso)))

# ========================================
# Step 3: Run All DE Methods
# ========================================
cat("\nStep 3: Running differential expression analysis...\n")
cat("  This may take a few minutes...\n")

results <- run_all_methods(
  dso,
  methods = c("DESeq2", "edgeR", "limma-voom")
)

cat(sprintf("  - Successfully ran %d methods\n", length(results)))

# Display summary for each method
for (method in names(results)) {
  n_genes <- nrow(results[[method]]$results)
  cat(sprintf("  - %s: %d genes analyzed\n", method, n_genes))
}

# ========================================
# Step 4: Filter Results
# ========================================
cat("\nStep 4: Filtering significant genes...\n")

padj_threshold <- 0.05
lfc_threshold <- 1

filtered_results <- lapply(results, function(res) {
  filter_results(
    res$results,
    padj_threshold = padj_threshold,
    lfc_threshold = lfc_threshold
  )
})

for (method in names(filtered_results)) {
  n_sig <- nrow(filtered_results[[method]])
  cat(sprintf("  - %s: %d significant genes (padj < %.2f, |log2FC| > %.1f)\n", 
              method, n_sig, padj_threshold, lfc_threshold))
}

# ========================================
# Step 5: Compare Methods
# ========================================
cat("\nStep 5: Comparing methods...\n")

comparison <- compare_methods(
  results,
  padj_threshold = padj_threshold,
  lfc_threshold = lfc_threshold
)

cat("\nSummary Statistics:\n")
print(comparison$summary_stats)

cat("\nOverlap Matrix:\n")
print(comparison$comparison_matrix)

cat(sprintf("\nConsensus genes (found by all methods): %d\n", 
            length(comparison$consensus_genes)))

# ========================================
# Step 6: Create Consensus Gene List
# ========================================
cat("\nStep 6: Creating consensus gene list...\n")

consensus <- create_consensus_list(
  results,
  min_methods = 2,  # Genes found by at least 2 methods
  padj_threshold = padj_threshold,
  lfc_threshold = lfc_threshold,
  aggregate_lfc = "median"
)

cat(sprintf("  - Found %d consensus genes (supported by >= 2 methods)\n", 
            nrow(consensus)))

cat("\nTop 10 Consensus Genes:\n")
print(head(consensus, 10))

# ========================================
# Step 7: Enrichment Analysis
# ========================================
cat("\nStep 7: Running enrichment analysis...\n")
cat("  Note: This requires internet connection and may take a few minutes\n")

# Only run enrichment if we have enough genes
if (nrow(consensus) > 10) {
  tryCatch({
    enrichment <- run_enrichment(
      gene_list = consensus$gene,
      gene_id_type = "ENSEMBL",
      organism = "human",
      analysis_type = "both",
      ontology = "BP",
      pvalue_cutoff = 0.05
    )
    
    if (!is.null(enrichment$GO)) {
      cat(sprintf("  - Found %d significant GO terms\n", enrichment$GO$n_significant))
      cat("\nTop 5 GO Terms:\n")
      print(head(enrichment$GO$summary[, c("Description", "pvalue", "Count")], 5))
    }
    
    if (!is.null(enrichment$KEGG)) {
      cat(sprintf("  - Found %d significant KEGG pathways\n", enrichment$KEGG$n_significant))
    }
  }, error = function(e) {
    cat("  - Enrichment analysis skipped (requires additional setup)\n")
    cat(sprintf("    Error: %s\n", e$message))
  })
}

# ========================================
# Step 8: Export Results
# ========================================
cat("\nStep 8: Exporting results...\n")

# Create output directory
output_dir <- "desandbox_output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Export all results
export_results(
  results,
  output_prefix = file.path(output_dir, "airway_analysis"),
  formats = c("csv", "tsv", "json", "rds")
)

# Export consensus genes
write.csv(
  consensus,
  file = file.path(output_dir, "consensus_genes.csv"),
  row.names = FALSE
)

cat(sprintf("  - Results exported to '%s/'\n", output_dir))

# ========================================
# Step 9: Cache Analysis
# ========================================
cat("\nStep 9: Caching results for reproducibility...\n")

cache_dir <- file.path(output_dir, "cache")
if (!dir.exists(cache_dir)) {
  dir.create(cache_dir, recursive = TRUE)
}

cache_path <- cache_analysis(
  results,
  cache_dir = cache_dir,
  cache_name = "airway_analysis",
  metadata = list(
    dataset = "airway",
    organism = "human",
    comparison = "dexamethasone vs control",
    date = Sys.Date(),
    n_genes = nrow(dso),
    n_samples = ncol(dso),
    padj_threshold = padj_threshold,
    lfc_threshold = lfc_threshold
  )
)

cat(sprintf("  - Results cached to: %s\n", cache_path))

# ========================================
# Summary
# ========================================
cat("\n=== Analysis Complete! ===\n\n")
cat("Summary:\n")
cat(sprintf("  - Dataset: airway (dexamethasone treatment)\n"))
cat(sprintf("  - Genes analyzed: %d\n", nrow(counts)))
cat(sprintf("  - Samples: %d (%d control, %d treated)\n", 
            ncol(counts),
            sum(metadata$condition == "untrt"),
            sum(metadata$condition == "trt")))
cat(sprintf("  - Methods: %s\n", paste(names(results), collapse = ", ")))
cat(sprintf("  - Consensus genes: %d\n", nrow(consensus)))
cat(sprintf("  - Output directory: %s/\n", output_dir))

cat("\nOutput files:\n")
cat("  - airway_analysis_DESeq2.csv\n")
cat("  - airway_analysis_edgeR.csv\n")
cat("  - airway_analysis_limma_voom.csv\n")
cat("  - airway_analysis_all.json\n")
cat("  - consensus_genes.csv\n")
cat("  - cache/airway_analysis.rds\n")

cat("\nTo reload cached results:\n")
cat("  cached_results <- load_cached_analysis('", cache_path, "')\n", sep = "")

cat("\n=== Done! ===\n")
