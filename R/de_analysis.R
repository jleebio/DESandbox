#' Run DESeq2 Analysis
#'
#' Performs differential expression analysis using DESeq2.
#'
#' @param desandbox_object A SummarizedExperiment object from create_desandbox_object()
#' @param design_formula A formula for the design matrix (default: ~condition)
#' @param covariates Character vector of covariate column names to include in design (default: NULL)
#' @param contrast A character vector specifying the contrast (default: NULL, uses first two levels)
#' @param alpha Significance level for FDR (default: 0.05)
#' @param lfc_threshold Log2 fold change threshold (default: 0)
#' @param parallel Whether to use parallel processing (default: FALSE)
#' @param ... Additional arguments passed to DESeq2::DESeq()
#'
#' @return A list containing DESeq2 results and metadata
#' @export
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom SummarizedExperiment assay colData
#' @importFrom BiocParallel bpparam
#'
#' @examples
#' # Create example data
#' counts <- matrix(rpois(1000, 10), nrow=100, ncol=10)
#' colnames(counts) <- paste0("sample", 1:10)
#' metadata <- data.frame(
#'   sample = colnames(counts),
#'   condition = factor(rep(c("control", "treatment"), each=5))
#' )
#' rownames(metadata) <- metadata$sample
#' 
#' # Create DESandbox object
#' dso <- create_desandbox_object(counts, metadata, "condition")
#' 
#' # Run DESeq2
#' results <- run_deseq2(dso)
#' head(results$results)
run_deseq2 <- function(desandbox_object,
                       design_formula = NULL,
                       covariates = NULL,
                       contrast = NULL,
                       alpha = 0.05,
                       lfc_threshold = 0,
                       parallel = FALSE,
                       ...) {
  
  # Extract data from SummarizedExperiment
  counts <- SummarizedExperiment::assay(desandbox_object, "counts")
  coldata <- as.data.frame(SummarizedExperiment::colData(desandbox_object))
  condition_col <- S4Vectors::metadata(desandbox_object)$condition_column
  
  # Set up design formula
  if (is.null(design_formula)) {
    if (!is.null(covariates)) {
      # Include covariates in design
      formula_terms <- c(covariates, condition_col)
      design_formula <- as.formula(paste0("~", paste(formula_terms, collapse = " + ")))
      message(sprintf("Using design formula: %s", deparse(design_formula)))
    } else {
      design_formula <- as.formula(paste0("~", condition_col))
    }
  }
  
  # Create DESeq2 object
  message("Creating DESeqDataSet...")
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = design_formula
  )
  
  # Run DESeq2
  message("Running DESeq2 analysis...")
  if (parallel) {
    dds <- DESeq2::DESeq(dds, parallel = TRUE, BPPARAM = BiocParallel::bpparam(), ...)
  } else {
    dds <- DESeq2::DESeq(dds, ...)
  }
  
  # Extract results
  message("Extracting results...")
  if (is.null(contrast)) {
    # Use default contrast (first two levels)
    res <- DESeq2::results(dds, alpha = alpha, lfcThreshold = lfc_threshold)
  } else {
    res <- DESeq2::results(dds, contrast = contrast, alpha = alpha, lfcThreshold = lfc_threshold)
  }
  
  # Convert to data.frame
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  # Reorder columns
  res_df <- res_df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  
  # Return structured results
  return(list(
    method = "DESeq2",
    results = res_df,
    dds_object = dds,
    res_object = res,
    parameters = list(
      design_formula = design_formula,
      covariates = covariates,
      contrast = contrast,
      alpha = alpha,
      lfc_threshold = lfc_threshold
    ),
    analysis_date = Sys.time()
  ))
}


#' Run edgeR Analysis
#'
#' Performs differential expression analysis using edgeR with quasi-likelihood F-test.
#'
#' @param desandbox_object A SummarizedExperiment object from create_desandbox_object()
#' @param design_formula A formula for the design matrix (default: ~condition)
#' @param covariates Character vector of covariate column names to include in design (default: NULL)
#' @param contrast A numeric contrast vector or column index (default: NULL, uses coef=2)
#' @param normalization Normalization method: "TMM", "RLE", "upperquartile", "none" (default: "TMM")
#' @param robust Whether to use robust estimation (default: TRUE)
#' @param ... Additional arguments passed to edgeR functions
#'
#' @return A list containing edgeR results and metadata
#' @export
#'
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmQLFit glmQLFTest topTags
#' @importFrom SummarizedExperiment assay colData
#' @importFrom stats model.matrix
#'
#' @examples
#' # See vignette for complete example
run_edger <- function(desandbox_object,
                      design_formula = NULL,
                      covariates = NULL,
                      contrast = NULL,
                      normalization = "TMM",
                      robust = TRUE,
                      ...) {
  
  # Extract data from SummarizedExperiment
  counts <- SummarizedExperiment::assay(desandbox_object, "counts")
  coldata <- as.data.frame(SummarizedExperiment::colData(desandbox_object))
  condition_col <- S4Vectors::metadata(desandbox_object)$condition_column
  
  # Set up design formula
  if (is.null(design_formula)) {
    if (!is.null(covariates)) {
      # Include covariates in design
      formula_terms <- c(covariates, condition_col)
      design_formula <- as.formula(paste0("~", paste(formula_terms, collapse = " + ")))
      message(sprintf("Using design formula: %s", deparse(design_formula)))
    } else {
      design_formula <- as.formula(paste0("~", condition_col))
    }
  }
  
  # Create design matrix
  design <- stats::model.matrix(design_formula, data = coldata)
  
  # Create DGEList object
  message("Creating DGEList...")
  dge <- edgeR::DGEList(counts = counts)
  
  # Normalize
  message(sprintf("Normalizing with %s method...", normalization))
  dge <- edgeR::calcNormFactors(dge, method = normalization)
  
  # Estimate dispersion
  message("Estimating dispersion...")
  dge <- edgeR::estimateDisp(dge, design, robust = robust)
  
  # Fit GLM
  message("Fitting quasi-likelihood GLM...")
  fit <- edgeR::glmQLFit(dge, design, robust = robust)
  
  # Perform test
  message("Performing QL F-test...")
  if (is.null(contrast)) {
    # Use default coefficient (usually the second column)
    qlf <- edgeR::glmQLFTest(fit, coef = 2)
  } else {
    qlf <- edgeR::glmQLFTest(fit, contrast = contrast)
  }
  
  # Extract top tags (all genes)
  tt <- edgeR::topTags(qlf, n = Inf, sort.by = "PValue")
  res_df <- as.data.frame(tt)
  
  # Standardize column names
  colnames(res_df)[colnames(res_df) == "logFC"] <- "log2FoldChange"
  colnames(res_df)[colnames(res_df) == "logCPM"] <- "baseMean"
  colnames(res_df)[colnames(res_df) == "PValue"] <- "pvalue"
  colnames(res_df)[colnames(res_df) == "FDR"] <- "padj"
  
  # Add gene column
  res_df$gene <- rownames(res_df)
  
  # Reorder columns to match standard format
  standard_cols <- c("gene", "baseMean", "log2FoldChange", "pvalue", "padj")
  other_cols <- setdiff(colnames(res_df), standard_cols)
  res_df <- res_df[, c(standard_cols, other_cols)]
  
  # Return structured results
  return(list(
    method = "edgeR",
    results = res_df,
    dge_object = dge,
    fit_object = fit,
    qlf_object = qlf,
    parameters = list(
      design_formula = design_formula,
      covariates = covariates,
      contrast = contrast,
      normalization = normalization,
      robust = robust
    ),
    analysis_date = Sys.time()
  ))
}


#' Run limma-voom Analysis
#'
#' Performs differential expression analysis using limma-voom.
#'
#' @param desandbox_object A SummarizedExperiment object from create_desandbox_object()
#' @param design_formula A formula for the design matrix (default: ~condition)
#' @param covariates Character vector of covariate column names to include in design (default: NULL)
#' @param contrast A numeric contrast vector or NULL (default: NULL, uses coef=2)
#' @param normalization Normalization method for voom: "TMM", "RLE", "upperquartile", "none" (default: "TMM")
#' @param plot Whether to generate voom mean-variance plot (default: FALSE)
#' @param robust Whether to use robust empirical Bayes (default: FALSE)
#' @param ... Additional arguments passed to limma functions
#'
#' @return A list containing limma-voom results and metadata
#' @export
#'
#' @importFrom limma voom lmFit eBayes topTable
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom SummarizedExperiment assay colData
#' @importFrom stats model.matrix
#'
#' @examples
#' # See vignette for complete example
run_limma_voom <- function(desandbox_object,
                           design_formula = NULL,
                           covariates = NULL,
                           contrast = NULL,
                           normalization = "TMM",
                           plot = FALSE,
                           robust = FALSE,
                           ...) {
  
  # Extract data from SummarizedExperiment
  counts <- SummarizedExperiment::assay(desandbox_object, "counts")
  coldata <- as.data.frame(SummarizedExperiment::colData(desandbox_object))
  condition_col <- S4Vectors::metadata(desandbox_object)$condition_column
  
  # Set up design formula
  if (is.null(design_formula)) {
    if (!is.null(covariates)) {
      # Include covariates in design
      formula_terms <- c(covariates, condition_col)
      design_formula <- as.formula(paste0("~", paste(formula_terms, collapse = " + ")))
      message(sprintf("Using design formula: %s", deparse(design_formula)))
    } else {
      design_formula <- as.formula(paste0("~", condition_col))
    }
  }
  
  # Create design matrix
  design <- stats::model.matrix(design_formula, data = coldata)
  
  # Create DGEList and normalize
  message("Creating DGEList and normalizing...")
  dge <- edgeR::DGEList(counts = counts)
  dge <- edgeR::calcNormFactors(dge, method = normalization)
  
  # Apply voom transformation
  message("Applying voom transformation...")
  v <- limma::voom(dge, design, plot = plot)
  
  # Fit linear model
  message("Fitting linear model...")
  fit <- limma::lmFit(v, design)
  
  # Apply empirical Bayes
  message("Applying empirical Bayes moderation...")
  fit <- limma::eBayes(fit, robust = robust)
  
  # Extract results
  message("Extracting results...")
  if (is.null(contrast)) {
    # Use default coefficient
    tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "P")
  } else {
    tt <- limma::topTable(fit, contrast = contrast, number = Inf, sort.by = "P")
  }
  
  # Standardize column names
  res_df <- as.data.frame(tt)
  colnames(res_df)[colnames(res_df) == "logFC"] <- "log2FoldChange"
  colnames(res_df)[colnames(res_df) == "AveExpr"] <- "baseMean"
  colnames(res_df)[colnames(res_df) == "P.Value"] <- "pvalue"
  colnames(res_df)[colnames(res_df) == "adj.P.Val"] <- "padj"
  
  # Add gene column
  res_df$gene <- rownames(res_df)
  
  # Reorder columns
  standard_cols <- c("gene", "baseMean", "log2FoldChange", "pvalue", "padj")
  other_cols <- setdiff(colnames(res_df), standard_cols)
  res_df <- res_df[, c(standard_cols, other_cols)]
  
  # Return structured results
  return(list(
    method = "limma-voom",
    results = res_df,
    dge_object = dge,
    voom_object = v,
    fit_object = fit,
    parameters = list(
      design_formula = design_formula,
      covariates = covariates,
      contrast = contrast,
      normalization = normalization,
      robust = robust
    ),
    analysis_date = Sys.time()
  ))
}


#' Run All DE Methods
#'
#' Convenience function to run DESeq2, edgeR, and limma-voom in parallel.
#'
#' @param desandbox_object A SummarizedExperiment object from create_desandbox_object()
#' @param methods Vector of methods to run: "DESeq2", "edgeR", "limma-voom" (default: all)
#' @param ... Additional arguments passed to individual method functions
#'
#' @return A list containing results from all methods
#' @export
#'
#' @examples
#' # Create example data
#' counts <- matrix(rpois(1000, 10), nrow=100, ncol=10)
#' colnames(counts) <- paste0("sample", 1:10)
#' metadata <- data.frame(
#'   sample = colnames(counts),
#'   condition = factor(rep(c("control", "treatment"), each=5))
#' )
#' rownames(metadata) <- metadata$sample
#' 
#' # Create DESandbox object and run all methods
#' dso <- create_desandbox_object(counts, metadata, "condition")
#' results <- run_all_methods(dso)
#' names(results)
run_all_methods <- function(desandbox_object,
                           methods = c("DESeq2", "edgeR", "limma-voom"),
                           ...) {
  
  results <- list()
  
  if ("DESeq2" %in% methods) {
    message("\n=== Running DESeq2 ===")
    results$DESeq2 <- tryCatch({
      run_deseq2(desandbox_object, ...)
    }, error = function(e) {
      warning(sprintf("DESeq2 failed: %s", e$message))
      NULL
    })
  }
  
  if ("edgeR" %in% methods) {
    message("\n=== Running edgeR ===")
    results$edgeR <- tryCatch({
      run_edger(desandbox_object, ...)
    }, error = function(e) {
      warning(sprintf("edgeR failed: %s", e$message))
      NULL
    })
  }
  
  if ("limma-voom" %in% methods) {
    message("\n=== Running limma-voom ===")
    results$`limma-voom` <- tryCatch({
      run_limma_voom(desandbox_object, ...)
    }, error = function(e) {
      warning(sprintf("limma-voom failed: %s", e$message))
      NULL
    })
  }
  
  # Remove NULL entries (failed methods)
  results <- results[!sapply(results, is.null)]
  
  if (length(results) == 0) {
    stop("All methods failed. Check your data and parameters.")
  }
  
  message(sprintf("\n=== Completed %d/%d methods successfully ===", 
                  length(results), length(methods)))
  
  return(results)
}
