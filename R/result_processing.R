#' Standardize DE Results
#'
#' Converts results from different DE methods into a common format.
#'
#' @param de_result A result object from run_deseq2(), run_edger(), or run_limma_voom()
#' @param method Method name (usually auto-detected from de_result)
#'
#' @return A standardized data.frame with consistent column names
#' @export
#'
#' @examples
#' # See vignette for complete example
standardize_results <- function(de_result, method = NULL) {
  
  if (is.null(method)) {
    method <- de_result$method
  }
  
  if (is.null(method)) {
    stop("Could not determine method. Please specify method argument.")
  }
  
  res_df <- de_result$results
  
  # Ensure required columns exist
  required_cols <- c("gene", "log2FoldChange", "pvalue", "padj")
  missing_cols <- setdiff(required_cols, colnames(res_df))
  
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }
  
  # Create standardized result
  standard_res <- data.frame(
    gene = res_df$gene,
    baseMean = if ("baseMean" %in% colnames(res_df)) res_df$baseMean else NA_real_,
    log2FoldChange = res_df$log2FoldChange,
    lfcSE = if ("lfcSE" %in% colnames(res_df)) res_df$lfcSE else NA_real_,
    stat = if ("stat" %in% colnames(res_df)) res_df$stat else 
           if ("t" %in% colnames(res_df)) res_df$t else 
           if ("F" %in% colnames(res_df)) res_df$F else NA_real_,
    pvalue = res_df$pvalue,
    padj = res_df$padj,
    method = method,
    stringsAsFactors = FALSE
  )
  
  # Remove rows with NA gene names
  standard_res <- standard_res[!is.na(standard_res$gene), ]
  
  # Sort by adjusted p-value
  standard_res <- standard_res[order(standard_res$padj, decreasing = FALSE), ]
  
  return(standard_res)
}


#' Filter DE Results
#'
#' Filters differential expression results based on thresholds.
#'
#' @param results A data.frame of DE results (standardized or raw)
#' @param padj_threshold Adjusted p-value threshold (default: 0.05)
#' @param lfc_threshold Absolute log2 fold change threshold (default: 1)
#' @param baseMean_threshold Minimum base mean expression (default: NULL)
#' @param direction Filter by direction: "both", "up", "down" (default: "both")
#'
#' @return A filtered data.frame
#' @export
#'
#' @examples
#' # See vignette for complete example
filter_results <- function(results,
                          padj_threshold = 0.05,
                          lfc_threshold = 1,
                          baseMean_threshold = NULL,
                          direction = "both") {
  
  # Check required columns
  if (!"padj" %in% colnames(results)) {
    stop("Results must contain 'padj' column")
  }
  
  if (!"log2FoldChange" %in% colnames(results)) {
    stop("Results must contain 'log2FoldChange' column")
  }
  
  # Apply filters
  filtered <- results
  
  # Filter by adjusted p-value
  filtered <- filtered[!is.na(filtered$padj) & filtered$padj < padj_threshold, ]
  
  # Filter by log2 fold change
  if (direction == "both") {
    filtered <- filtered[!is.na(filtered$log2FoldChange) & 
                          abs(filtered$log2FoldChange) >= lfc_threshold, ]
  } else if (direction == "up") {
    filtered <- filtered[!is.na(filtered$log2FoldChange) & 
                          filtered$log2FoldChange >= lfc_threshold, ]
  } else if (direction == "down") {
    filtered <- filtered[!is.na(filtered$log2FoldChange) & 
                          filtered$log2FoldChange <= -lfc_threshold, ]
  } else {
    stop("direction must be 'both', 'up', or 'down'")
  }
  
  # Filter by base mean if specified
  if (!is.null(baseMean_threshold) && "baseMean" %in% colnames(filtered)) {
    filtered <- filtered[!is.na(filtered$baseMean) & 
                          filtered$baseMean >= baseMean_threshold, ]
  }
  
  message(sprintf("Filtered to %d significant genes (from %d total)", 
                  nrow(filtered), nrow(results)))
  
  return(filtered)
}


#' Export Results to Multiple Formats
#'
#' Exports DE results to CSV, TSV, or JSON formats for downstream use.
#'
#' @param results A data.frame or list of results
#' @param output_prefix Output file prefix (without extension)
#' @param formats Vector of formats: "csv", "tsv", "json", "rds" (default: all)
#' @param include_metadata Whether to include analysis metadata in JSON (default: TRUE)
#'
#' @return A list of output file paths
#' @export
#'
#' @importFrom utils write.table
#' @importFrom jsonlite toJSON
#'
#' @examples
#' # See vignette for complete example
export_results <- function(results,
                          output_prefix,
                          formats = c("csv", "tsv", "json", "rds"),
                          include_metadata = TRUE) {
  
  output_files <- character()
  
  # Handle list of results (multiple methods)
  if (is.list(results) && !is.data.frame(results)) {
    
    for (format in formats) {
      if (format == "json") {
        # Export all results as single JSON
        output_file <- paste0(output_prefix, "_all.", format)
        
        # Prepare data for JSON export
        export_data <- lapply(names(results), function(method) {
          res <- results[[method]]
          
          if (include_metadata) {
            # Convert formula to string for JSON serialization
            params <- res$parameters
            if (!is.null(params$design_formula) && inherits(params$design_formula, "formula")) {
              params$design_formula <- deparse(params$design_formula)
            }
            
            list(
              method = method,
              results = res$results,
              parameters = params,
              analysis_date = as.character(res$analysis_date)
            )
          } else {
            res$results
          }
        })
        names(export_data) <- names(results)
        
        json_str <- jsonlite::toJSON(export_data, pretty = TRUE, auto_unbox = TRUE)
        writeLines(json_str, output_file)
        output_files <- c(output_files, output_file)
        message(sprintf("Exported to %s", output_file))
        
      } else if (format == "rds") {
        # Save complete R object
        output_file <- paste0(output_prefix, "_all.", format)
        saveRDS(results, output_file)
        output_files <- c(output_files, output_file)
        message(sprintf("Saved R object to %s", output_file))
        
      } else {
        # Export each method separately for CSV/TSV
        for (method in names(results)) {
          res_df <- results[[method]]$results
          safe_method <- gsub("-", "_", method)
          output_file <- paste0(output_prefix, "_", safe_method, ".", format)
          
          if (format == "csv") {
            utils::write.table(res_df, output_file, sep = ",", 
                             row.names = FALSE, quote = FALSE)
          } else if (format == "tsv") {
            utils::write.table(res_df, output_file, sep = "\t", 
                             row.names = FALSE, quote = FALSE)
          }
          
          output_files <- c(output_files, output_file)
          message(sprintf("Exported %s results to %s", method, output_file))
        }
      }
    }
    
  } else {
    # Handle single result data.frame
    res_df <- if (is.data.frame(results)) results else results$results
    
    for (format in formats) {
      output_file <- paste0(output_prefix, ".", format)
      
      if (format == "csv") {
        utils::write.table(res_df, output_file, sep = ",", 
                         row.names = FALSE, quote = FALSE)
      } else if (format == "tsv") {
        utils::write.table(res_df, output_file, sep = "\t", 
                         row.names = FALSE, quote = FALSE)
      } else if (format == "json") {
        export_data <- if (is.data.frame(results)) {
          list(results = res_df)
        } else if (include_metadata) {
          # Convert formula to string for JSON serialization
          params <- results$parameters
          if (!is.null(params$design_formula) && inherits(params$design_formula, "formula")) {
            params$design_formula <- deparse(params$design_formula)
          }
          
          list(
            method = results$method,
            results = res_df,
            parameters = params,
            analysis_date = as.character(results$analysis_date)
          )
        } else {
          res_df
        }
        
        json_str <- jsonlite::toJSON(export_data, pretty = TRUE, auto_unbox = TRUE)
        writeLines(json_str, output_file)
      } else if (format == "rds") {
        saveRDS(results, output_file)
      }
      
      output_files <- c(output_files, output_file)
      message(sprintf("Exported to %s", output_file))
    }
  }
  
  return(invisible(output_files))
}
