# DESandbox Package Architecture

## Overview

DESandbox is designed as a modular, extensible R package following Bioconductor best practices. The architecture emphasizes:

- **Modularity**: Clear separation of concerns across functional modules
- **Testability**: Each module has well-defined inputs/outputs
- **Reproducibility**: Built-in caching and parameter tracking
- **API Integration**: JSON export and standardized data structures

## Directory Structure

```
DESandbox/
├── DESCRIPTION           # Package metadata and dependencies
├── NAMESPACE            # Exported functions and imports
├── LICENSE              # MIT license
├── README.md            # User-facing documentation
│
├── R/                   # Core package code
│   ├── data_validation.R      # Input validation and object creation
│   ├── de_analysis.R          # DE method wrappers (DESeq2, edgeR, limma)
│   ├── result_processing.R    # Standardization, filtering, export
│   ├── comparison.R           # Cross-method comparison
│   ├── enrichment.R           # GO/KEGG enrichment analysis
│   └── caching.R             # Result caching and reproducibility
│
├── man/                 # Auto-generated documentation (roxygen2)
│
├── tests/               # Unit tests
│   ├── testthat.R
│   └── testthat/
│       └── test-desandbox.R
│
├── vignettes/          # User guides and tutorials
│   └── getting_started.Rmd
│
└── inst/               # Additional package data
    └── extdata/        # Example datasets
```

## Module Descriptions

### 1. Data Validation Module (`data_validation.R`)

**Purpose**: Ensure data quality before analysis

**Key Functions**:
- `validate_counts()`: Check count matrix format, filter low-count genes
- `validate_metadata()`: Verify sample metadata consistency
- `create_desandbox_object()`: Create SummarizedExperiment object

**Design Principles**:
- Fail fast with clear error messages
- Provide warnings for sub-optimal data
- Return structured validation results

**Data Flow**:
```
Raw counts + metadata 
  → validate_counts() 
  → validate_metadata() 
  → create_desandbox_object() 
  → SummarizedExperiment
```

### 2. DE Analysis Module (`de_analysis.R`)

**Purpose**: Unified interface to three DE methods

**Key Functions**:
- `run_deseq2()`: Wrapper for DESeq2 analysis
- `run_edger()`: Wrapper for edgeR quasi-likelihood F-test
- `run_limma_voom()`: Wrapper for limma-voom analysis
- `run_all_methods()`: Execute all methods in parallel

**Design Principles**:
- Consistent API across all methods
- Error handling to prevent pipeline failures
- Return standardized result structure

**Result Structure**:
```r
list(
  method = "DESeq2",
  results = data.frame(...),  # Standardized columns
  <method>_object = ...,       # Original object for advanced users
  parameters = list(...),      # Analysis parameters
  analysis_date = ...          # Timestamp
)
```

### 3. Result Processing Module (`result_processing.R`)

**Purpose**: Standardize and export results

**Key Functions**:
- `standardize_results()`: Convert to common format
- `filter_results()`: Apply significance thresholds
- `export_results()`: Export to CSV/TSV/JSON/RDS

**Standardized Result Schema**:
```r
data.frame(
  gene = character,           # Gene identifier
  baseMean = numeric,         # Mean expression
  log2FoldChange = numeric,   # Effect size
  lfcSE = numeric,           # LFC standard error
  stat = numeric,            # Test statistic
  pvalue = numeric,          # Raw p-value
  padj = numeric,            # Adjusted p-value (FDR)
  method = character         # Analysis method
)
```

### 4. Comparison Module (`comparison.R`)

**Purpose**: Compare results across methods

**Key Functions**:
- `compare_methods()`: Comprehensive cross-method comparison
- `compute_overlap()`: Detailed overlap analysis
- `create_consensus_list()`: Generate consensus gene list

**Comparison Metrics**:
- Overlap counts (intersection)
- Jaccard similarity index
- Method-specific genes (unique findings)
- Consensus genes (all methods agree)

### 5. Enrichment Module (`enrichment.R`)

**Purpose**: Gene set enrichment analysis

**Key Functions**:
- `run_enrichment()`: GO and KEGG enrichment
- `run_enrichment_from_de()`: Convenience wrapper for DE results
- `compare_enrichment()`: Cross-method pathway comparison

**Supported Analyses**:
- GO (Gene Ontology): BP, MF, CC
- KEGG pathways
- Custom gene sets (future)

### 6. Caching Module (`caching.R`)

**Purpose**: Reproducibility and performance

**Key Functions**:
- `cache_analysis()`: Save results with metadata
- `load_cached_analysis()`: Load cached results
- `list_cache_files()`: Browse cache
- `clear_cache()`: Clean old cache files

**Cache Structure**:
```r
list(
  results = <actual_results>,
  cache_metadata = list(
    cached_date = ...,
    desandbox_version = ...,
    r_version = ...,
    user_metadata = list(...)
  )
)
```

## Data Structures

### Primary Input: SummarizedExperiment

```r
SummarizedExperiment(
  assays = list(counts = matrix),
  colData = DataFrame(metadata),
  metadata = list(
    condition_column = "condition",
    creation_date = ...,
    desandbox_version = ...
  )
)
```

### Primary Output: Standardized Results List

```r
list(
  method = character,
  results = data.frame,
  <method>_object = S4 object,
  parameters = list,
  analysis_date = POSIXct
)
```

## Function Naming Conventions

- `validate_*()`: Data validation functions
- `run_*()`: Analysis execution functions
- `create_*()`: Object creation functions
- `compute_*()`: Computation/calculation functions
- `filter_*()`: Data filtering functions
- `export_*()`: Data export functions

## Error Handling Strategy

1. **Validation errors**: Stop with clear message
2. **Analysis errors**: Catch and return NULL with warning
3. **Export errors**: Catch and report failed exports
4. **Enrichment errors**: Warn but don't stop pipeline

## Testing Strategy

- **Unit tests**: Each function tested independently
- **Integration tests**: End-to-end workflow tests
- **Mock data**: Simulated datasets for reproducibility
- **Edge cases**: Empty data, single samples, etc.

## Dependencies

### Core Dependencies
- **DESeq2**: DESeq2 differential expression
- **edgeR**: edgeR differential expression
- **limma**: limma-voom differential expression
- **SummarizedExperiment**: Data containers
- **clusterProfiler**: Enrichment analysis

### Utility Dependencies
- **digest**: Hashing for caching
- **jsonlite**: JSON export
- **BiocParallel**: Parallel processing

## Performance Considerations

1. **Memory**: Use sparse matrices for large datasets (future)
2. **Speed**: Parallel processing for DESeq2
3. **Caching**: Avoid re-running expensive analyses
4. **Filtering**: Pre-filter low-count genes early

## Future Enhancements

1. **Additional methods**: Add more DE tools (sleuth, ballgown)
2. **Visualization**: Built-in plotting functions
3. **Batch effects**: Integration with sva/RUVSeq
4. **Single-cell**: Extend to scRNA-seq (Seurat, Scanpy compatibility)
5. **Custom gene sets**: User-defined enrichment databases
6. **Report generation**: Automated HTML/PDF reports

## Best Practices for Contributors

1. **Documentation**: All functions must have roxygen2 docs
2. **Testing**: Add tests for new functionality
3. **Code style**: Follow Bioconductor style guide
4. **Backward compatibility**: Deprecate, don't break
5. **Performance**: Profile before optimizing
6. **Dependencies**: Minimize new dependencies
