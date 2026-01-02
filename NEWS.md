# DESandbox 0.99.0

## New Features

* Initial Bioconductor submission
* Unified interface for differential expression analysis using DESeq2, edgeR, and limma-voom
* Standardized output format across all methods for easy comparison
* Publication-quality visualization tools:
  - Volcano plots with customizable thresholds and gene labeling
  - Venn diagrams for 2-way and 3-way method comparisons
  - Bar plots for DE gene count comparisons
  - Nature Genetics publication styling with threshold annotations
* Cross-method comparison utilities:
  - Compute overlap and consensus between methods
  - Identify method-specific and shared DE genes
  - Generate comparison statistics
* Gene set enrichment analysis integration (GO and KEGG)
* Covariate support in design formulas for all methods
* Result caching system for reproducibility
* Export functionality to multiple formats (CSV, TSV, JSON, RDS)
* Comprehensive vignette with real data examples (airway dataset)
* Full test coverage with 25 unit tests

## Documentation

* All 17 exported functions fully documented with roxygen2
* 34 help files with examples
* Complete vignette demonstrating full workflow
* Additional guides for covariates and visualizations
