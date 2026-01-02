#' Compare DE Results Across Methods
#'
#' Compares differential expression results from multiple methods.
#'
#' @param results_list A list of DE results from different methods
#' @param padj_threshold Adjusted p-value threshold for significance (default: 0.05)
#' @param lfc_threshold Absolute log2 fold change threshold (default: 1)
#'
#' @return A list containing comparison statistics and data
#' @export
#'
#' @examples
#' # Create mock results from multiple methods
#' results_list <- list(
#'   DESeq2 = list(results = data.frame(
#'     gene = paste0("Gene", 1:100),
#'     log2FoldChange = rnorm(100, 0, 2),
#'     padj = runif(100)
#'   )),
#'   edgeR = list(results = data.frame(
#'     gene = paste0("Gene", 1:100),
#'     log2FoldChange = rnorm(100, 0, 2),
#'     padj = runif(100)
#'   ))
#' )
#' 
#' # Compare methods
#' comparison <- compare_methods(results_list)
#' print(comparison$summary)
compare_methods <- function(results_list,
                           padj_threshold = 0.05,
                           lfc_threshold = 1) {
  
  if (!is.list(results_list) || length(results_list) < 2) {
    stop("results_list must be a list with at least 2 methods")
  }
  
  method_names <- names(results_list)
  if (is.null(method_names)) {
    method_names <- paste0("Method", seq_along(results_list))
    names(results_list) <- method_names
  }
  
  # Extract and filter significant genes for each method
  sig_genes_list <- list()
  all_results <- list()
  
  for (method in method_names) {
    res <- results_list[[method]]
    res_df <- if (is.data.frame(res)) res else res$results
    
    # Filter for significant genes
    sig <- filter_results(
      res_df,
      padj_threshold = padj_threshold,
      lfc_threshold = lfc_threshold,
      direction = "both"
    )
    
    sig_genes_list[[method]] <- sig$gene
    all_results[[method]] <- res_df
  }
  
  # Compute overlap statistics
  n_methods <- length(method_names)
  comparison_matrix <- matrix(0, nrow = n_methods, ncol = n_methods,
                             dimnames = list(method_names, method_names))
  
  for (i in 1:n_methods) {
    for (j in 1:n_methods) {
      if (i <= j) {
        overlap <- length(intersect(sig_genes_list[[i]], sig_genes_list[[j]]))
        comparison_matrix[i, j] <- overlap
        comparison_matrix[j, i] <- overlap
      }
    }
  }
  
  # Calculate Jaccard index for each pair
  jaccard_matrix <- matrix(0, nrow = n_methods, ncol = n_methods,
                          dimnames = list(method_names, method_names))
  
  for (i in 1:n_methods) {
    for (j in 1:n_methods) {
      if (i < j) {
        intersection <- length(intersect(sig_genes_list[[i]], sig_genes_list[[j]]))
        union <- length(union(sig_genes_list[[i]], sig_genes_list[[j]]))
        jaccard <- if (union > 0) intersection / union else 0
        jaccard_matrix[i, j] <- jaccard
        jaccard_matrix[j, i] <- jaccard
      } else if (i == j) {
        jaccard_matrix[i, j] <- 1.0
      }
    }
  }
  
  # Count genes found by different numbers of methods
  all_genes <- unique(unlist(sig_genes_list))
  genes_by_method_count <- sapply(all_genes, function(gene) {
    sum(sapply(sig_genes_list, function(x) gene %in% x))
  })
  
  method_count_table <- table(genes_by_method_count)
  
  # Create summary statistics
  summary_stats <- data.frame(
    method = method_names,
    n_significant = sapply(sig_genes_list, length),
    n_upregulated = sapply(method_names, function(m) {
      res_df <- all_results[[m]]
      sig_genes <- sig_genes_list[[m]]
      sum(res_df$gene %in% sig_genes & res_df$log2FoldChange > 0, na.rm = TRUE)
    }),
    n_downregulated = sapply(method_names, function(m) {
      res_df <- all_results[[m]]
      sig_genes <- sig_genes_list[[m]]
      sum(res_df$gene %in% sig_genes & res_df$log2FoldChange < 0, na.rm = TRUE)
    }),
    stringsAsFactors = FALSE
  )
  
  # Find consensus genes (found by all methods)
  consensus_genes <- Reduce(intersect, sig_genes_list)
  
  # Find method-specific genes (unique to each method)
  method_specific <- lapply(method_names, function(m) {
    others <- setdiff(method_names, m)
    other_genes <- unique(unlist(sig_genes_list[others]))
    setdiff(sig_genes_list[[m]], other_genes)
  })
  names(method_specific) <- method_names
  
  return(list(
    comparison_matrix = comparison_matrix,
    jaccard_matrix = jaccard_matrix,
    summary_stats = summary_stats,
    method_count_table = method_count_table,
    consensus_genes = consensus_genes,
    method_specific_genes = method_specific,
    sig_genes_list = sig_genes_list,
    parameters = list(
      padj_threshold = padj_threshold,
      lfc_threshold = lfc_threshold
    )
  ))
}


#' Compute Gene Overlap Between Methods
#'
#' Computes detailed overlap statistics and Venn diagram data.
#'
#' @param results_list A list of DE results from different methods
#' @param padj_threshold Adjusted p-value threshold (default: 0.05)
#' @param lfc_threshold Absolute log2 fold change threshold (default: 1)
#' @param direction Filter by direction: "both", "up", "down" (default: "both")
#'
#' @return A list containing overlap data and statistics
#' @export
#'
#' @examples
#' # Create mock results
#' results_list <- list(
#'   DESeq2 = list(results = data.frame(
#'     gene = paste0("Gene", 1:100),
#'     log2FoldChange = rnorm(100, 0, 2),
#'     padj = runif(100)
#'   )),
#'   edgeR = list(results = data.frame(
#'     gene = paste0("Gene", 1:100),
#'     log2FoldChange = rnorm(100, 0, 2),
#'     padj = runif(100)
#'   ))
#' )
#' 
#' # Compute overlap
#' overlap <- compute_overlap(results_list)
#' print(overlap$summary)
compute_overlap <- function(results_list,
                           padj_threshold = 0.05,
                           lfc_threshold = 1,
                           direction = "both") {
  
  if (!is.list(results_list)) {
    stop("results_list must be a list of DE results")
  }
  
  method_names <- names(results_list)
  if (is.null(method_names)) {
    method_names <- paste0("Method", seq_along(results_list))
    names(results_list) <- method_names
  }
  
  # Extract significant genes for each method
  sig_genes_list <- list()
  
  for (method in method_names) {
    res <- results_list[[method]]
    res_df <- if (is.data.frame(res)) res else res$results
    
    # Filter for significant genes
    sig <- filter_results(
      res_df,
      padj_threshold = padj_threshold,
      lfc_threshold = lfc_threshold,
      direction = direction
    )
    
    sig_genes_list[[method]] <- sig$gene
  }
  
  # Calculate all pairwise overlaps
  n_methods <- length(method_names)
  pairwise_overlaps <- list()
  
  for (i in 1:(n_methods - 1)) {
    for (j in (i + 1):n_methods) {
      method1 <- method_names[i]
      method2 <- method_names[j]
      
      overlap_genes <- intersect(sig_genes_list[[method1]], sig_genes_list[[method2]])
      unique_to_1 <- setdiff(sig_genes_list[[method1]], sig_genes_list[[method2]])
      unique_to_2 <- setdiff(sig_genes_list[[method2]], sig_genes_list[[method1]])
      
      pair_name <- paste(method1, method2, sep = "_vs_")
      pairwise_overlaps[[pair_name]] <- list(
        method1 = method1,
        method2 = method2,
        n_overlap = length(overlap_genes),
        n_unique_method1 = length(unique_to_1),
        n_unique_method2 = length(unique_to_2),
        overlap_genes = overlap_genes,
        unique_to_method1 = unique_to_1,
        unique_to_method2 = unique_to_2,
        jaccard_index = length(overlap_genes) / length(union(sig_genes_list[[method1]], 
                                                              sig_genes_list[[method2]]))
      )
    }
  }
  
  # Calculate overlap for all methods combined
  if (n_methods >= 2) {
    all_genes <- unique(unlist(sig_genes_list))
    
    overlap_counts <- sapply(all_genes, function(gene) {
      sum(sapply(sig_genes_list, function(x) gene %in% x))
    })
    
    overlap_by_count <- lapply(1:n_methods, function(k) {
      all_genes[overlap_counts == k]
    })
    names(overlap_by_count) <- paste0("found_in_", 1:n_methods, "_methods")
  } else {
    overlap_by_count <- NULL
  }
  
  # Create a summary data.frame for easy viewing
  overlap_summary <- do.call(rbind, lapply(names(pairwise_overlaps), function(pair_name) {
    pair <- pairwise_overlaps[[pair_name]]
    data.frame(
      comparison = pair_name,
      method1 = pair$method1,
      method2 = pair$method2,
      n_overlap = pair$n_overlap,
      n_unique_method1 = pair$n_unique_method1,
      n_unique_method2 = pair$n_unique_method2,
      jaccard_index = round(pair$jaccard_index, 3),
      stringsAsFactors = FALSE
    )
  }))
  
  return(list(
    pairwise_overlaps = pairwise_overlaps,
    overlap_by_count = overlap_by_count,
    overlap_summary = overlap_summary,
    sig_genes_list = sig_genes_list,
    parameters = list(
      padj_threshold = padj_threshold,
      lfc_threshold = lfc_threshold,
      direction = direction
    )
  ))
}


#' Create Consensus Gene List
#'
#' Creates a consensus list of DE genes based on agreement across methods.
#'
#' @param results_list A list of DE results from different methods
#' @param min_methods Minimum number of methods that must agree (default: 2)
#' @param padj_threshold Adjusted p-value threshold (default: 0.05)
#' @param lfc_threshold Absolute log2 fold change threshold (default: 1)
#' @param aggregate_lfc How to aggregate log2FC: "mean", "median", "min", "max" (default: "median")
#'
#' @return A data.frame with consensus DE genes
#' @export
#'
#' @examples
#' # Create mock results from 3 methods
#' results_list <- list(
#'   DESeq2 = list(results = data.frame(
#'     gene = paste0("Gene", 1:100),
#'     log2FoldChange = rnorm(100, 0, 2),
#'     padj = runif(100)
#'   )),
#'   edgeR = list(results = data.frame(
#'     gene = paste0("Gene", 1:100),
#'     log2FoldChange = rnorm(100, 0, 2),
#'     padj = runif(100)
#'   )),
#'   limma = list(results = data.frame(
#'     gene = paste0("Gene", 1:100),
#'     log2FoldChange = rnorm(100, 0, 2),
#'     padj = runif(100)
#'   ))
#' )
#' 
#' # Create consensus requiring 2+ methods
#' consensus <- create_consensus_list(results_list, min_methods = 2)
#' head(consensus)
create_consensus_list <- function(results_list,
                                  min_methods = 2,
                                  padj_threshold = 0.05,
                                  lfc_threshold = 1,
                                  aggregate_lfc = "median") {
  
  method_names <- names(results_list)
  
  # Get significant genes from each method
  sig_genes_list <- list()
  all_results_list <- list()
  
  for (method in method_names) {
    res <- results_list[[method]]
    res_df <- if (is.data.frame(res)) res else res$results
    all_results_list[[method]] <- res_df
    
    sig <- filter_results(
      res_df,
      padj_threshold = padj_threshold,
      lfc_threshold = lfc_threshold
    )
    
    sig_genes_list[[method]] <- sig$gene
  }
  
  # Find all genes and count support
  all_genes <- unique(unlist(sig_genes_list))
  
  consensus_data <- lapply(all_genes, function(gene) {
    # Count how many methods found this gene
    n_support <- sum(sapply(sig_genes_list, function(x) gene %in% x))
    
    # Extract stats from each method
    lfc_values <- c()
    padj_values <- c()
    methods_found <- c()
    
    for (method in method_names) {
      res_df <- all_results_list[[method]]
      gene_row <- res_df[res_df$gene == gene, ]
      
      if (nrow(gene_row) > 0 && gene %in% sig_genes_list[[method]]) {
        lfc_values <- c(lfc_values, gene_row$log2FoldChange[1])
        padj_values <- c(padj_values, gene_row$padj[1])
        methods_found <- c(methods_found, method)
      }
    }
    
    # Aggregate statistics
    if (aggregate_lfc == "mean") {
      agg_lfc <- mean(lfc_values, na.rm = TRUE)
    } else if (aggregate_lfc == "median") {
      agg_lfc <- median(lfc_values, na.rm = TRUE)
    } else if (aggregate_lfc == "min") {
      agg_lfc <- min(abs(lfc_values), na.rm = TRUE) * sign(median(lfc_values, na.rm = TRUE))
    } else if (aggregate_lfc == "max") {
      agg_lfc <- max(abs(lfc_values), na.rm = TRUE) * sign(median(lfc_values, na.rm = TRUE))
    } else {
      agg_lfc <- median(lfc_values, na.rm = TRUE)
    }
    
    data.frame(
      gene = gene,
      n_methods_support = n_support,
      log2FoldChange = agg_lfc,
      mean_padj = mean(padj_values, na.rm = TRUE),
      min_padj = min(padj_values, na.rm = TRUE),
      methods = paste(methods_found, collapse = ","),
      stringsAsFactors = FALSE
    )
  })
  
  consensus_df <- do.call(rbind, consensus_data)
  
  # Filter by minimum method support
  consensus_df <- consensus_df[consensus_df$n_methods_support >= min_methods, ]
  
  # Sort by support and then by p-value
  consensus_df <- consensus_df[order(-consensus_df$n_methods_support, 
                                     consensus_df$min_padj), ]
  
  message(sprintf("Found %d consensus genes (supported by >= %d methods)", 
                  nrow(consensus_df), min_methods))
  
  return(consensus_df)
}
