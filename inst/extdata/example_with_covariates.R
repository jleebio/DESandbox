#!/usr/bin/env Rscript

# Example: Using Covariates in DESandbox
# This demonstrates how to include batch effects or other covariates
# in your differential expression analysis

# Load required packages
library(airway)
library(SummarizedExperiment)

# Source DESandbox functions (for development testing)
source("R/data_validation.R")
source("R/de_analysis.R")
source("R/result_processing.R")
source("R/comparison.R")

# Load the airway dataset
data(airway)
counts_matrix <- assay(airway, "counts")
metadata_df <- as.data.frame(colData(airway))

# The airway dataset has a "cell" column (cell line)
# We can treat this as a covariate/batch effect to control for
print("Sample metadata:")
print(metadata_df[, c("dex", "cell", "SampleName")])

# Create DESandbox object
dso <- create_desandbox_object(
  counts = counts_matrix,
  metadata = metadata_df,
  condition_column = "dex"
)

# ==========================================
# Example 1: Analysis WITHOUT covariates
# ==========================================
cat("\n=== Example 1: Standard analysis (no covariates) ===\n")
cat("Design formula: ~dex\n\n")

result_standard <- run_deseq2(dso)
sig_genes_standard <- nrow(filter_results(
  result_standard$results,
  padj_threshold = 0.05,
  lfc_threshold = 1
))
cat(sprintf("Significant genes (no covariates): %d\n", sig_genes_standard))

# ==========================================
# Example 2: Analysis WITH cell line covariate
# ==========================================
cat("\n=== Example 2: Analysis controlling for cell line ===\n")
cat("Design formula: ~cell + dex\n\n")

result_with_covariate <- run_deseq2(
  dso,
  covariates = "cell"  # Control for cell line differences
)
sig_genes_covariate <- nrow(filter_results(
  result_with_covariate$results,
  padj_threshold = 0.05,
  lfc_threshold = 1
))
cat(sprintf("Significant genes (with cell covariate): %d\n", sig_genes_covariate))

# ==========================================
# Example 3: Multiple covariates
# ==========================================
# If you had multiple covariates (e.g., batch, sex, age), you can include them all:
# result_multi <- run_deseq2(dso, covariates = c("batch", "sex", "age"))
# This creates design: ~batch + sex + age + condition

cat("\n=== Example 3: Running all methods with covariates ===\n")
all_results <- run_all_methods(dso, covariates = "cell")

# Compare results
cat("\nSignificant genes by method (with cell covariate):\n")
for (method_name in names(all_results)) {
  if (!is.null(all_results[[method_name]])) {
    sig_count <- nrow(filter_results(
      all_results[[method_name]]$results,
      padj_threshold = 0.05,
      lfc_threshold = 1
    ))
    cat(sprintf("  %s: %d genes\n", method_name, sig_count))
  }
}

# ==========================================
# Example 4: Paired samples design
# ==========================================
# For paired samples (e.g., before/after treatment on same subjects):
# If you have a "subject_id" column identifying paired samples:
# 
# paired_result <- run_deseq2(
#   dso,
#   covariates = "subject_id"  # Accounts for subject-specific baseline
# )
# Design: ~subject_id + condition

cat("\n=== Summary ===\n")
cat("Covariates help account for known sources of variation:\n")
cat("  - Batch effects (processing date, sequencing run)\n")
cat("  - Biological covariates (cell line, sex, age)\n")
cat("  - Paired samples (subject ID, time points)\n")
cat("  - Technical replicates (library prep batch)\n\n")
cat("Design formula construction:\n")
cat("  - No covariates: ~condition\n")
cat("  - Single covariate: ~batch + condition\n")
cat("  - Multiple covariates: ~batch + sex + age + condition\n")
cat("  - Always put the condition of interest LAST\n")
