#' Validate Count Matrix
#'
#' Validates that the count matrix meets requirements for DE analysis.
#'
#' @param counts A numeric matrix or data.frame with genes in rows and samples in columns
#' @param min_counts Minimum total counts per gene (default: 10)
#' @param min_samples Minimum number of samples a gene must be expressed in (default: 3)
#' @param allow_negative Whether to allow negative values (default: FALSE)
#'
#' @return A list with validation status and cleaned count matrix
#' @export
#'
#' @examples
#' counts <- matrix(rpois(1000, 10), nrow=100, ncol=10)
#' result <- validate_counts(counts)
validate_counts <- function(counts, 
                           min_counts = 10, 
                           min_samples = 3,
                           allow_negative = FALSE) {
  
  # Check if counts is a matrix or can be coerced to one
  if (!is.matrix(counts) && !is.data.frame(counts)) {
    stop("counts must be a matrix or data.frame")
  }
  
  # Convert to matrix if data.frame
  if (is.data.frame(counts)) {
    gene_names <- rownames(counts)
    counts <- as.matrix(counts)
    rownames(counts) <- gene_names
  }
  
  # Check for numeric values
  if (!is.numeric(counts)) {
    stop("counts must contain numeric values")
  }
  
  # Check for negative values
  if (!allow_negative && any(counts < 0, na.rm = TRUE)) {
    stop("counts contains negative values. Set allow_negative=TRUE to override.")
  }
  
  # Check for NA values
  if (any(is.na(counts))) {
    warning("counts contains NA values. These will be replaced with 0.")
    counts[is.na(counts)] <- 0
  }
  
  # Check for non-integer values (warning only)
  if (any(counts != floor(counts), na.rm = TRUE)) {
    warning("counts contains non-integer values. Consider rounding.")
  }
  
  # Filter low-count genes
  gene_sums <- rowSums(counts)
  expressed_samples <- rowSums(counts > 0)
  
  keep <- gene_sums >= min_counts & expressed_samples >= min_samples
  
  n_filtered <- sum(!keep)
  if (n_filtered > 0) {
    message(sprintf("Filtered %d genes with low counts", n_filtered))
  }
  
  filtered_counts <- counts[keep, , drop = FALSE]
  
  # Check dimensions
  if (nrow(filtered_counts) < 100) {
    warning("Very few genes remaining after filtering. Consider relaxing filters.")
  }
  
  if (ncol(filtered_counts) < 3) {
    stop("At least 3 samples required for DE analysis")
  }
  
  return(list(
    valid = TRUE,
    counts = filtered_counts,
    n_genes = nrow(filtered_counts),
    n_samples = ncol(filtered_counts),
    n_filtered = n_filtered,
    filter_params = list(
      min_counts = min_counts,
      min_samples = min_samples
    )
  ))
}


#' Validate Sample Metadata
#'
#' Validates that sample metadata is properly formatted and matches the count matrix.
#'
#' @param metadata A data.frame with sample information
#' @param counts A count matrix (for matching samples)
#' @param condition_column Name of the column containing condition/group information
#' @param required_columns Vector of required column names
#'
#' @return A list with validation status and validated metadata
#' @export
#'
#' @examples
#' metadata <- data.frame(
#'   sample = paste0("sample", 1:10),
#'   condition = rep(c("control", "treatment"), each=5)
#' )
#' result <- validate_metadata(metadata, condition_column = "condition")
validate_metadata <- function(metadata, 
                              counts = NULL,
                              condition_column = "condition",
                              required_columns = NULL) {
  
  # Check if metadata is a data.frame
  if (!is.data.frame(metadata)) {
    stop("metadata must be a data.frame")
  }
  
  # Check if required columns exist
  if (!is.null(required_columns)) {
    missing_cols <- setdiff(required_columns, colnames(metadata))
    if (length(missing_cols) > 0) {
      stop(sprintf("Missing required columns: %s", 
                   paste(missing_cols, collapse = ", ")))
    }
  }
  
  # Check if condition column exists
  if (!condition_column %in% colnames(metadata)) {
    stop(sprintf("Condition column '%s' not found in metadata", condition_column))
  }
  
  # Convert condition to factor if not already
  if (!is.factor(metadata[[condition_column]])) {
    metadata[[condition_column]] <- factor(metadata[[condition_column]])
    message(sprintf("Converted '%s' to factor", condition_column))
  }
  
  # Check for at least 2 conditions
  n_conditions <- nlevels(metadata[[condition_column]])
  if (n_conditions < 2) {
    stop("At least 2 conditions required for DE analysis")
  }
  
  # Check for balanced design (warning only)
  condition_counts <- table(metadata[[condition_column]])
  if (any(condition_counts < 2)) {
    warning("Some conditions have fewer than 2 replicates. Results may be unreliable.")
  }
  
  # Match with count matrix if provided
  if (!is.null(counts)) {
    # Ensure rownames exist
    if (is.null(rownames(metadata))) {
      if ("sample" %in% colnames(metadata)) {
        rownames(metadata) <- metadata$sample
      } else {
        stop("metadata must have rownames or a 'sample' column")
      }
    }
    
    # Check sample matching
    count_samples <- colnames(counts)
    meta_samples <- rownames(metadata)
    
    if (!all(count_samples %in% meta_samples)) {
      missing <- setdiff(count_samples, meta_samples)
      stop(sprintf("Count samples not in metadata: %s", 
                   paste(missing, collapse = ", ")))
    }
    
    if (!all(meta_samples %in% count_samples)) {
      extra <- setdiff(meta_samples, count_samples)
      warning(sprintf("Metadata samples not in counts: %s", 
                      paste(extra, collapse = ", ")))
    }
    
    # Reorder metadata to match counts
    metadata <- metadata[count_samples, , drop = FALSE]
  }
  
  return(list(
    valid = TRUE,
    metadata = metadata,
    n_samples = nrow(metadata),
    n_conditions = n_conditions,
    condition_column = condition_column,
    conditions = levels(metadata[[condition_column]])
  ))
}


#' Create DESandbox Object
#'
#' Creates a standardized S4 object containing counts and metadata for DE analysis.
#'
#' @param counts A numeric matrix with genes in rows and samples in columns
#' @param metadata A data.frame with sample information
#' @param condition_column Name of the column containing condition information
#' @param validate Whether to validate inputs (default: TRUE)
#' @param ... Additional arguments passed to validation functions
#'
#' @return A SummarizedExperiment object
#' @export
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#'
#' @examples
#' counts <- matrix(rpois(1000, 10), nrow=100, ncol=10)
#' colnames(counts) <- paste0("sample", 1:10)
#' metadata <- data.frame(
#'   sample = colnames(counts),
#'   condition = rep(c("control", "treatment"), each=5)
#' )
#' rownames(metadata) <- metadata$sample
#' dso <- create_desandbox_object(counts, metadata, "condition")
create_desandbox_object <- function(counts, 
                                   metadata, 
                                   condition_column = "condition",
                                   validate = TRUE,
                                   ...) {
  
  if (validate) {
    # Validate counts
    count_validation <- validate_counts(counts, ...)
    counts <- count_validation$counts
    
    # Validate metadata
    meta_validation <- validate_metadata(metadata, counts, condition_column, ...)
    metadata <- meta_validation$metadata
  }
  
  # Create SummarizedExperiment object
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(metadata)
  )
  
  # Add metadata about creation
  S4Vectors::metadata(se)$condition_column <- condition_column
  S4Vectors::metadata(se)$creation_date <- Sys.time()
  S4Vectors::metadata(se)$desandbox_version <- tryCatch(
    as.character(utils::packageVersion("DESandbox")),
    error = function(e) "0.1.0-dev"
  )
  
  return(se)
}
