# DESandbox Quick Reference

## Installation

```r
BiocManager::install(c("DESeq2", "edgeR", "limma", "clusterProfiler", "org.Hs.eg.db"))
# devtools::install_github("yourusername/DESandbox")
```

## Basic Workflow

```r
library(DESandbox)

# 1. Load data
counts <- read.csv("counts.csv", row.names = 1)
metadata <- read.csv("metadata.csv", row.names = 1)

# 2. Create object
dso <- create_desandbox_object(counts, metadata, "condition")

# 3. Run analysis
results <- run_all_methods(dso)

# 4. Compare methods
comparison <- compare_methods(results, padj_threshold = 0.05, lfc_threshold = 1)

# 5. Get consensus genes
consensus <- create_consensus_list(results, min_methods = 2)

# 6. Run enrichment
enrichment <- compare_enrichment(results, organism = "human")

# 7. Export
export_results(results, "output", formats = c("csv", "json"))
```

## Function Quick Reference

### Data Validation
```r
validate_counts(counts, min_counts = 10, min_samples = 3)
validate_metadata(metadata, condition_column = "condition")
create_desandbox_object(counts, metadata, "condition")
```

### Run Analysis
```r
# All methods
run_all_methods(dso, methods = c("DESeq2", "edgeR", "limma-voom"))

# Individual methods
run_deseq2(dso, alpha = 0.05, lfc_threshold = 0)
run_edger(dso, normalization = "TMM", robust = TRUE)
run_limma_voom(dso, normalization = "TMM", plot = FALSE)

# With covariates (batch correction, paired samples)
run_deseq2(dso, covariates = "batch")                      # Single covariate
run_edger(dso, covariates = c("batch", "sex"))            # Multiple covariates
run_limma_voom(dso, covariates = "subject_id")            # Paired design
run_all_methods(dso, covariates = c("batch", "cell"))     # All methods with covariates
```

### Filter & Process
```r
filter_results(results, padj_threshold = 0.05, lfc_threshold = 1, direction = "both")
standardize_results(de_result)
```

### Comparison
```r
compare_methods(results_list, padj_threshold = 0.05, lfc_threshold = 1)
compute_overlap(results_list, padj_threshold = 0.05, lfc_threshold = 1)
create_consensus_list(results_list, min_methods = 2, aggregate_lfc = "median")
```

### Enrichment
```r
run_enrichment(gene_list, gene_id_type = "SYMBOL", organism = "human", 
               analysis_type = "both", ontology = "BP")
run_enrichment_from_de(de_results, padj_threshold = 0.05, lfc_threshold = 1)
compare_enrichment(results_list, organism = "human")
```

### Export & Caching
```r
export_results(results, "prefix", formats = c("csv", "tsv", "json", "rds"))
cache_analysis(results, cache_dir = "./cache")
load_cached_analysis(cache_path)
list_cache_files("./cache")
clear_cache(cache_dir = "./cache", older_than = 30)
```

## Common Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `padj_threshold` | Adjusted p-value cutoff | 0.05 |
| `lfc_threshold` | Absolute log2 fold change cutoff | 1 |
| `direction` | "both", "up", or "down" | "both" |
| `min_counts` | Minimum total counts per gene | 10 |
| `min_samples` | Minimum samples a gene must be in | 3 |
| `organism` | "human", "mouse", or "rat" | "human" |
| `analysis_type` | "GO", "KEGG", or "both" | "both" |
| `ontology` | "BP", "MF", "CC", or "ALL" | "BP" |

## Result Structure

```r
# DE result structure
list(
  method = "DESeq2",
  results = data.frame(
    gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
  ),
  dds_object = ...,  # Original DESeq2 object
  parameters = list(...),
  analysis_date = ...
)

# Comparison structure
list(
  comparison_matrix = matrix,  # Overlap counts
  jaccard_matrix = matrix,     # Jaccard similarity
  summary_stats = data.frame,  # Per-method statistics
  consensus_genes = character, # Genes found by all methods
  method_specific_genes = list # Unique genes per method
)

# Enrichment structure
list(
  GO = list(result, summary, n_significant),
  KEGG = list(result, summary, n_significant),
  parameters = list(...)
)
```

## Typical Thresholds

### Conservative
```r
padj_threshold = 0.01
lfc_threshold = 2
min_methods = 3  # All methods must agree
```

### Standard
```r
padj_threshold = 0.05
lfc_threshold = 1
min_methods = 2  # At least 2 methods agree
```

### Exploratory
```r
padj_threshold = 0.1
lfc_threshold = 0.5
min_methods = 1  # Any method
```

## Troubleshooting

### Too few genes after filtering
```r
# Relax filtering thresholds
count_validation <- validate_counts(counts, min_counts = 5, min_samples = 2)

# Or adjust DE thresholds
filter_results(results, padj_threshold = 0.1, lfc_threshold = 0.5)
```

### Method fails
```r
# Run methods individually to isolate issue
deseq2_results <- tryCatch(
  run_deseq2(dso),
  error = function(e) {
    message(e$message)
    NULL
  }
)
```

### Enrichment returns no results
```r
# Check gene ID format
table(nchar(gene_list))  # ENTREZ IDs are usually 4-6 digits

# Try different ID type
enrichment <- run_enrichment(
  gene_list,
  gene_id_type = "ENSEMBL",  # Try ENSEMBL instead of SYMBOL
  organism = "human"
)

# Increase p-value threshold
enrichment <- run_enrichment(
  gene_list,
  pvalue_cutoff = 0.1,
  qvalue_cutoff = 0.3
)
```

## Example Datasets

### Airway dataset (Bioconductor)
```r
library(airway)
data(airway)
counts <- assay(airway)
metadata <- as.data.frame(colData(airway))
```

### Simulated data
```r
# Create test data
counts <- matrix(rnbinom(10000, mu = 100, size = 10), nrow = 1000, ncol = 10)
rownames(counts) <- paste0("Gene", 1:1000)
colnames(counts) <- paste0("Sample", 1:10)

metadata <- data.frame(
  sample = colnames(counts),
  condition = factor(rep(c("control", "treatment"), each = 5))
)
rownames(metadata) <- metadata$sample
```

## Best Practices

1. ✅ Always validate inputs before analysis
2. ✅ Run multiple methods and compare
3. ✅ Use consensus genes for downstream analysis
4. ✅ Cache large analyses
5. ✅ Export results in multiple formats
6. ✅ Document analysis parameters
7. ❌ Don't skip validation
8. ❌ Don't rely on a single method
9. ❌ Don't ignore filtering statistics
10. ❌ Don't forget to save parameters

## Getting Help

```r
# Function help
?create_desandbox_object
?run_all_methods
?compare_methods

# Package vignette
vignette("getting_started", package = "DESandbox")

# Check for errors
get_errors()  # If available in your environment
```

## Resources

- GitHub: https://github.com/yourusername/DESandbox
- Documentation: See `man/` directory
- Vignette: See `vignettes/getting_started.Rmd`
- Architecture: See `ARCHITECTURE.md`
- API Integration: See `API_INTEGRATION.md`
