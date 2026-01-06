# DESandbox

**Unified Differential Expression Analysis in R**

Compare DESeq2, edgeR, and limma-voom results with one simple interface.

## ✨ Why DESandbox?

- **Easy to use**: Same commands for all three DE methods
- **Compare results**: See how different methods agree
- **Beautiful plots**: Publication-ready visualizations
- **Find consensus**: Identify genes supported by multiple methods

## 📦 Installation

```r
# Install dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "edgeR", "limma"))

# Install DESandbox (coming soon on Bioconductor)
devtools::install_github("jleebio/DESandbox")
```

## 🚀 Getting Started

Quick start with the airway example (no vignette needed):

```r
# Install required packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("airway", "DESeq2", "edgeR", "limma"))

library(DESandbox)
library(airway)
library(SummarizedExperiment)

# Load example data
data("airway", package = "airway")

# Prepare counts and metadata
counts <- as.matrix(assay(airway))
metadata <- as.data.frame(colData(airway))
metadata$dex <- stats::relevel(factor(metadata$dex), ref = "untrt")

# Create DESandbox object
dso <- create_desandbox_object(counts, metadata, condition_column = "dex")

# Run DESeq2 (you can also run edgeR/limma-voom)
de_deseq2 <- run_deseq2(dso, alpha = 0.05, lfc_threshold = 1)

# Standardize and plot volcano
std_df <- standardize_results(de_deseq2)
std_df$gene_id <- std_df$gene
std_df$log2FC  <- std_df$log2FoldChange
plot_volcano(std_df, facet_by_method = FALSE, label_top_n = 10,
             title = "airway: DESeq2 Volcano")
```

## 📚 Learn More

See the **vignettes** folder for detailed guides:
- `getting_started.Rmd` - Complete workflow example
- Additional guides in `inst/doc/`

## 🎯 Quick Example (Run All Methods)

Prefer to compare methods in one go?

```r
results <- run_all_methods(dso)            # runs DESeq2, edgeR, limma-voom
comparison <- compare_methods(results)     # overlaps and consensus

# Example: plot DESeq2 volcano from combined results
de2_std <- standardize_results(results$DESeq2)
de2_std$gene_id <- de2_std$gene; de2_std$log2FC <- de2_std$log2FoldChange
plot_volcano(de2_std, facet_by_method = FALSE)
```

## 🔗 Links

- **Documentation**: See `vignettes/` folder (or use the quick start above)
- **Report Issues**: GitHub Issues
- **Bioconductor**: Coming soon

## 📄 License

MIT License - See LICENSE file

---

**Have questions?** Check the vignettes first - they cover everything!
