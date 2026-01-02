# DESandbox: Unified Differential Expression Analysis

A comprehensive R package providing a unified interface for running and comparing differential expression analyses using DESeq2, edgeR, and limma-voom.

## Features

- **Unified Interface**: Consistent API across three popular DE methods
- **Standardized Results**: Common output format for easy comparison
- **Method Comparison**: Built-in tools to compare and overlap results
- **Publication-Quality Visualizations**: Volcano plots, Venn diagrams, and comparison charts
- **Covariate Support**: Control for batch effects and paired sample designs
- **Gene Set Enrichment**: Integrated GO and KEGG pathway analysis
- **Reproducibility**: Caching system and parameter tracking
- **API-Ready**: JSON export and structured outputs for web integration

## Installation

```r
# Install from GitHub (development version)
# devtools::install_github("yourusername/DESandbox")

# Or install dependencies manually
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "DESeq2", "edgeR", "limma", 
    "SummarizedExperiment", "clusterProfiler",
    "org.Hs.eg.db"
))
```

## Quick Start

### Complete Working Example with Airway Dataset

```r
library(DESandbox)

# Install and load the airway dataset from Bioconductor
# BiocManager::install("airway")
library(airway)
data(airway)

# Extract counts and metadata
counts <- assay(airway)
metadata <- as.data.frame(colData(airway))

# Prepare metadata for analysis
metadata$condition <- metadata$dex
rownames(metadata) <- metadata$Run

# Filter low-count genes for faster example
counts <- counts[rowSums(counts) > 100, ]

# Create DESandbox object
dso <- create_desandbox_object(
    counts = counts,
    metadata = metadata,
    condition_column = "condition"
)

# Run all three methods (DESeq2, edgeR, limma-voom)
results <- run_all_methods(dso)

# Compare methods
comparison <- compare_methods(
    results,
    padj_threshold = 0.05,
    lfc_threshold = 1
)

# View summary statistics
print(comparison$summary_stats)
#>   method n_significant n_upregulated n_downregulated
#> 1 DESeq2           450           250             200
#> 2  edgeR           435           245             190
#> 3 limma-voom       442           248             194

# Get consensus genes (found by at least 2 methods)
consensus <- create_consensus_list(
    results,
    min_methods = 2
)

print(head(consensus))
#>         gene n_methods_support log2FoldChange  mean_padj   min_padj
#> 1 ENSG00000...                 3           2.45  0.000123  0.000089
#> 2 ENSG00000...                 3          -1.87  0.000456  0.000234

# Run enrichment analysis on consensus genes
enrichment <- run_enrichment(
    gene_list = consensus$gene,
    gene_id_type = "ENSEMBL",
    organism = "human",
    analysis_type = "both"
)

# View top GO terms
print(head(enrichment$GO$summary, 3))

# Export results in multiple formats
export_results(
    results,
    output_prefix = "airway_analysis",
    formats = c("csv", "json", "rds")
)

# Cache results for reproducibility
cache_path <- cache_analysis(
    results,
    cache_dir = "./cache",
    metadata = list(
        project = "Airway DEX Treatment",
        date = Sys.Date()
    )
)
```

### Quick Start with Your Own Data

```r
library(DESandbox)

# Load your data
counts <- read.csv("counts.csv", row.names = 1)
metadata <- read.csv("metadata.csv", row.names = 1)

# Run complete analysis in 3 lines
dso <- create_desandbox_object(counts, metadata, "condition")
results <- run_all_methods(dso)
consensus <- create_consensus_list(results, min_methods = 2)

# Export
export_results(results, "my_analysis", formats = c("csv", "json"))
```

### Controlling for Covariates (Batch Effects, Paired Samples)

DESandbox supports including covariates in your analysis to control for confounding factors:

```r
# Include covariates (e.g., batch effects, cell line, sex)
results_with_batch <- run_deseq2(
    dso,
    covariates = "batch"  # Design: ~batch + condition
)

# Multiple covariates
results_multi <- run_deseq2(
    dso,
    covariates = c("batch", "sex", "age")  # Design: ~batch + sex + age + condition
)

# Paired samples design (e.g., before/after on same subjects)
results_paired <- run_deseq2(
    dso,
    covariates = "subject_id"  # Accounts for subject baseline
)

# All methods support covariates
all_results <- run_all_methods(dso, covariates = "batch")
```

**When to use covariates:**
- Batch effects (sequencing run, processing date)
- Biological factors (cell line, sex, age, tissue type)
- Paired/matched samples (subject ID, time points)
- Technical replicates (library prep batch)

**Note:** Always put your condition of interest LAST in the design formula (automatically handled by DESandbox).

### Visualizations

DESandbox provides publication-quality visualizations:

```r
# Prepare combined results data frame
combined_df <- data.frame(rbind(
  data.frame(
    method = "DESeq2",
    gene_id = results$DESeq2$results$gene,
    log2FC = results$DESeq2$results$log2FoldChange,
    padj = results$DESeq2$results$padj
  ),
  data.frame(
    method = "edgeR",
    gene_id = results$edgeR$results$gene,
    log2FC = results$edgeR$results$log2FoldChange,
    padj = results$edgeR$results$padj
  )
))

# Volcano plots
p_volcano <- plot_volcano(
  combined_df,
  facet_by_method = TRUE,
  label_top_n = 10
)
print(p_volcano)

# Venn diagrams for gene overlap
p_venn_up <- plot_venn_de(combined_df, direction = "up")
p_venn_down <- plot_venn_de(combined_df, direction = "down")

# Get gene lists from Venn analysis
gene_sets <- plot_venn_de(combined_df, return_type = "sets")

# Comparison bar chart
p_counts <- plot_de_counts(combined_df)

# Save plots
ggplot2::ggsave("volcano.png", p_volcano, width = 12, height = 4, dpi = 300)
ggplot2::ggsave("venn_up.png", p_venn_up, width = 8, height = 8, dpi = 300)
```

See [VISUALIZATION_GUIDE.md](VISUALIZATION_GUIDE.md) for comprehensive examples.

### Run Complete Example Scripts

Full working examples are included in the package:

```r
# Basic analysis example
basic_example <- system.file(
    "extdata", 
    "example_airway_analysis.R", 
    package = "DESandbox"
)
source(basic_example)

# Covariate usage example
covariate_example <- system.file(
    "extdata",
    "example_with_covariates.R",
    package = "DESandbox"
)
source(covariate_example)

# Visualization example
viz_example <- system.file(
    "extdata",
    "example_visualizations.R",
    package = "DESandbox"
)
source(viz_example)
```

Or from the command line:

```bash
Rscript inst/extdata/example_airway_analysis.R
Rscript inst/extdata/example_visualizations.R
```

This will:
- Load the airway dataset
- Run all three DE methods
- Compare results and find consensus genes
- Run enrichment analysis
- Export results to `desandbox_output/`
- Cache results for reproducibility

## Typical Workflow

### 1. Data Preparation and Validation

```r
# Validate your count matrix
count_validation <- validate_counts(
    counts,
    min_counts = 10,
    min_samples = 3
)

# Validate metadata
meta_validation <- validate_metadata(
    metadata,
    counts = count_validation$counts,
    condition_column = "condition"
)
```

### 2. Create DESandbox Object

```r
dso <- create_desandbox_object(
    counts = counts,
    metadata = metadata,
    condition_column = "condition",
    validate = TRUE
)
```

### 3. Run Individual Methods

```r
# DESeq2
deseq2_results <- run_deseq2(
    dso,
    alpha = 0.05,
    lfc_threshold = 0
)

# edgeR
edger_results <- run_edger(
    dso,
    normalization = "TMM",
    robust = TRUE
)

# limma-voom
limma_results <- run_limma_voom(
    dso,
    normalization = "TMM",
    plot = FALSE
)
```

### 4. Filter and Compare Results

```r
# Filter results
sig_genes <- filter_results(
    deseq2_results$results,
    padj_threshold = 0.05,
    lfc_threshold = 1,
    direction = "both"
)

# Compare all methods
comparison <- compare_methods(
    list(
        DESeq2 = deseq2_results,
        edgeR = edger_results,
        limma = limma_results
    )
)

# View overlap statistics
print(comparison$summary_stats)
print(comparison$jaccard_matrix)
```

### 5. Enrichment Analysis

```r
# Run enrichment on consensus genes
enrichment <- run_enrichment(
    gene_list = consensus$gene,
    gene_id_type = "SYMBOL",
    organism = "human",
    analysis_type = "both",
    ontology = "BP"
)

# View top GO terms
head(enrichment$GO$summary)

# View top KEGG pathways
head(enrichment$KEGG$summary)
```

### 6. Export and Cache

```r
# Export results in multiple formats
export_results(
    results,
    output_prefix = "results/de_analysis",
    formats = c("csv", "tsv", "json", "rds")
)

# Cache for later use
cache_path <- cache_analysis(
    results,
    cache_dir = "./cache",
    metadata = list(
        project = "My Project",
        date = Sys.Date()
    )
)

# Load cached results later
cached_results <- load_cached_analysis(cache_path)
```

## Package Architecture

### Core Modules

- **data_validation.R**: Input validation and DESandbox object creation
- **de_analysis.R**: Wrappers for DESeq2, edgeR, and limma-voom
- **result_processing.R**: Standardization, filtering, and export
- **comparison.R**: Cross-method comparison and consensus analysis
- **enrichment.R**: GO and KEGG pathway enrichment
- **caching.R**: Result caching and reproducibility tracking

### Key Functions

| Function | Purpose |
|----------|---------|
| `create_desandbox_object()` | Create standardized input object |
| `run_all_methods()` | Run DESeq2, edgeR, and limma-voom |
| `compare_methods()` | Compare results across methods |
| `create_consensus_list()` | Generate consensus gene list |
| `run_enrichment()` | Perform GO/KEGG enrichment |
| `export_results()` | Export to CSV/JSON/RDS |
| `cache_analysis()` | Cache results for reproducibility |

## API Integration

DESandbox is designed for easy integration with web APIs (e.g., Django REST Framework):

```r
# Export results as JSON for API consumption
export_results(
    results,
    output_prefix = "api_output",
    formats = "json",
    include_metadata = TRUE
)

# The JSON structure is API-friendly:
# {
#   "DESeq2": {
#     "method": "DESeq2",
#     "results": [...],
#     "parameters": {...},
#     "analysis_date": "2026-01-02"
#   },
#   ...
# }
```

## Best Practices

1. **Always validate inputs** before analysis
2. **Use meaningful thresholds** (padj < 0.05, |log2FC| > 1)
3. **Compare multiple methods** to ensure robust results
4. **Cache large analyses** for reproducibility
5. **Document parameters** using the built-in metadata system
6. **Export in multiple formats** for different downstream tools

## Advanced Features

### Custom Contrasts

```r
# Specify custom contrasts
results <- run_deseq2(
    dso,
    contrast = c("condition", "treatment", "control")
)
```

### Parallel Processing

```r
# Enable parallel processing for DESeq2
results <- run_deseq2(dso, parallel = TRUE)
```

### Direction-Specific Analysis

```r
# Filter only upregulated genes
upregulated <- filter_results(
    results$DESeq2$results,
    direction = "up"
)
```

## Contributing

Contributions are welcome! Please follow Bioconductor coding standards.

## License

MIT License - see LICENSE file for details.

## Citation

If you use DESandbox in your research, please cite:

```
DESandbox: A unified interface for differential expression analysis
```

## Support

For issues and questions, please open a GitHub issue.
