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

**New to DESandbox?** Start with the vignette:

```r
library(DESandbox)
vignette("Getting Started with DESandbox")
```

This shows you:
- How to prepare your data
- How to run DE analysis
- How to compare methods
- How to create visualizations

## 📚 Learn More

See the **vignettes** folder for detailed guides:
- `getting_started.Rmd` - Complete workflow example
- Additional guides in `inst/doc/`

## 🎯 Quick Example

```r
library(DESandbox)

# Prepare your data
dso <- create_desandbox_object(
    counts = your_counts,
    metadata = your_metadata,
    condition_column = "condition"
)

# Run all three methods at once
results <- run_all_methods(dso)

# Compare the results
comparison <- compare_methods(results)

# Create a volcano plot
plot_volcano(results[[1]]$results)
```

## 🔗 Links

- **Documentation**: See vignettes (run `vignette("Getting Started with DESandbox")`)
- **Report Issues**: GitHub Issues
- **Bioconductor**: Coming soon

## 📄 License

MIT License - See LICENSE file

---

**Have questions?** Check the vignettes first - they cover everything!
