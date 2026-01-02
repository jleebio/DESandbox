#' Run Gene Set Enrichment Analysis
#'
#' Performs GO and/or KEGG enrichment analysis on DE gene lists.
#'
#' @param gene_list Character vector of gene IDs (or named vector with log2FC)
#' @param gene_id_type Type of gene IDs: "SYMBOL", "ENSEMBL", "ENTREZID" (default: "SYMBOL")
#' @param organism Organism database: "human", "mouse", "rat" (default: "human")
#' @param analysis_type Type of enrichment: "GO", "KEGG", "both" (default: "both")
#' @param ontology GO ontology: "BP", "MF", "CC", "ALL" (default: "BP")
#' @param pvalue_cutoff P-value cutoff (default: 0.05)
#' @param qvalue_cutoff Q-value cutoff (default: 0.2)
#' @param background_genes Optional background gene set (universe)
#' @param min_gs_size Minimum gene set size (default: 10)
#' @param max_gs_size Maximum gene set size (default: 500)
#'
#' @return A list containing enrichment results
#' @export
#'
#' @importFrom clusterProfiler enrichGO enrichKEGG
#' @importFrom AnnotationDbi mapIds
#'
#' @examples
#' # See vignette for complete example
run_enrichment <- function(gene_list,
                          gene_id_type = "SYMBOL",
                          organism = "human",
                          analysis_type = "both",
                          ontology = "BP",
                          pvalue_cutoff = 0.05,
                          qvalue_cutoff = 0.2,
                          background_genes = NULL,
                          min_gs_size = 10,
                          max_gs_size = 500) {
  
  # Select organism database
  if (organism == "human") {
    orgdb <- "org.Hs.eg.db"
    kegg_organism <- "hsa"
  } else if (organism == "mouse") {
    orgdb <- "org.Mm.eg.db"
    kegg_organism <- "mmu"
  } else if (organism == "rat") {
    orgdb <- "org.Rn.eg.db"
    kegg_organism <- "rno"
  } else {
    stop("Unsupported organism. Use 'human', 'mouse', or 'rat'")
  }
  
  # Load organism database
  if (!requireNamespace(orgdb, quietly = TRUE)) {
    stop(sprintf("Please install %s: BiocManager::install('%s')", orgdb, orgdb))
  }
  
  # Convert gene IDs to ENTREZ if needed
  if (gene_id_type != "ENTREZID") {
    message(sprintf("Converting %s to ENTREZID...", gene_id_type))
    
    gene_names <- if (is.null(names(gene_list))) gene_list else names(gene_list)
    
    entrez_ids <- AnnotationDbi::mapIds(
      get(orgdb),
      keys = gene_names,
      column = "ENTREZID",
      keytype = gene_id_type,
      multiVals = "first"
    )
    
    # Remove NAs
    entrez_ids <- entrez_ids[!is.na(entrez_ids)]
    
    if (length(entrez_ids) == 0) {
      stop("No genes could be mapped to ENTREZID. Check your gene IDs.")
    }
    
    message(sprintf("Mapped %d/%d genes to ENTREZID", 
                    length(entrez_ids), length(gene_list)))
    
    # If gene_list was named (with fold changes), preserve that
    if (!is.null(names(gene_list))) {
      gene_list <- gene_list[names(entrez_ids)]
      names(gene_list) <- as.character(entrez_ids)
    } else {
      gene_list <- as.character(entrez_ids)
    }
  }
  
  # Convert background genes if provided
  if (!is.null(background_genes)) {
    if (gene_id_type != "ENTREZID") {
      bg_entrez <- AnnotationDbi::mapIds(
        get(orgdb),
        keys = background_genes,
        column = "ENTREZID",
        keytype = gene_id_type,
        multiVals = "first"
      )
      background_genes <- as.character(bg_entrez[!is.na(bg_entrez)])
    }
  }
  
  results <- list()
  
  # GO Enrichment
  if (analysis_type %in% c("GO", "both")) {
    message("Running GO enrichment analysis...")
    
    go_result <- tryCatch({
      clusterProfiler::enrichGO(
        gene = if (is.null(names(gene_list))) gene_list else names(gene_list),
        OrgDb = orgdb,
        keyType = "ENTREZID",
        ont = ontology,
        pAdjustMethod = "BH",
        pvalueCutoff = pvalue_cutoff,
        qvalueCutoff = qvalue_cutoff,
        universe = background_genes,
        minGSSize = min_gs_size,
        maxGSSize = max_gs_size,
        readable = TRUE
      )
    }, error = function(e) {
      warning(sprintf("GO enrichment failed: %s", e$message))
      NULL
    })
    
    if (!is.null(go_result)) {
      results$GO <- list(
        result = go_result,
        summary = as.data.frame(go_result),
        n_significant = nrow(as.data.frame(go_result))
      )
      message(sprintf("Found %d significant GO terms", results$GO$n_significant))
    }
  }
  
  # KEGG Enrichment
  if (analysis_type %in% c("KEGG", "both")) {
    message("Running KEGG pathway enrichment...")
    
    kegg_result <- tryCatch({
      clusterProfiler::enrichKEGG(
        gene = if (is.null(names(gene_list))) gene_list else names(gene_list),
        organism = kegg_organism,
        keyType = "ncbi-geneid",
        pAdjustMethod = "BH",
        pvalueCutoff = pvalue_cutoff,
        qvalueCutoff = qvalue_cutoff,
        universe = background_genes,
        minGSSize = min_gs_size,
        maxGSSize = max_gs_size
      )
    }, error = function(e) {
      warning(sprintf("KEGG enrichment failed: %s", e$message))
      NULL
    })
    
    if (!is.null(kegg_result)) {
      results$KEGG <- list(
        result = kegg_result,
        summary = as.data.frame(kegg_result),
        n_significant = nrow(as.data.frame(kegg_result))
      )
      message(sprintf("Found %d significant KEGG pathways", results$KEGG$n_significant))
    }
  }
  
  # Add metadata
  results$parameters <- list(
    gene_id_type = gene_id_type,
    organism = organism,
    analysis_type = analysis_type,
    ontology = ontology,
    pvalue_cutoff = pvalue_cutoff,
    qvalue_cutoff = qvalue_cutoff,
    n_genes = length(gene_list),
    analysis_date = Sys.time()
  )
  
  return(results)
}


#' Run Enrichment on DE Results
#'
#' Convenience wrapper to run enrichment directly on filtered DE results.
#'
#' @param de_results DE results data.frame or result object
#' @param padj_threshold Adjusted p-value threshold (default: 0.05)
#' @param lfc_threshold Absolute log2 fold change threshold (default: 1)
#' @param direction Direction filter: "both", "up", "down" (default: "both")
#' @param use_lfc Whether to use log2FC values for ranked analysis (default: FALSE)
#' @param ... Additional arguments passed to run_enrichment()
#'
#' @return Enrichment results
#' @export
#'
#' @examples
#' # See vignette for complete example
run_enrichment_from_de <- function(de_results,
                                   padj_threshold = 0.05,
                                   lfc_threshold = 1,
                                   direction = "both",
                                   use_lfc = FALSE,
                                   ...) {
  
  # Extract data.frame if needed
  res_df <- if (is.data.frame(de_results)) de_results else de_results$results
  
  # Filter for significant genes
  sig_genes <- filter_results(
    res_df,
    padj_threshold = padj_threshold,
    lfc_threshold = lfc_threshold,
    direction = direction
  )
  
  if (nrow(sig_genes) == 0) {
    stop("No significant genes found with current thresholds")
  }
  
  message(sprintf("Running enrichment on %d significant genes", nrow(sig_genes)))
  
  # Prepare gene list
  if (use_lfc) {
    # Create named vector with log2FC values
    gene_list <- sig_genes$log2FoldChange
    names(gene_list) <- sig_genes$gene
    gene_list <- sort(gene_list, decreasing = TRUE)
  } else {
    # Just use gene names
    gene_list <- sig_genes$gene
  }
  
  # Optionally use all genes as background
  background <- res_df$gene
  
  # Run enrichment
  enrichment_results <- run_enrichment(
    gene_list = gene_list,
    background_genes = background,
    ...
  )
  
  return(enrichment_results)
}


#' Compare Enrichment Across Methods
#'
#' Runs enrichment for each method and compares results.
#'
#' @param results_list List of DE results from different methods
#' @param padj_threshold Adjusted p-value threshold (default: 0.05)
#' @param lfc_threshold Absolute log2 fold change threshold (default: 1)
#' @param ... Additional arguments passed to run_enrichment()
#'
#' @return List of enrichment results for each method plus comparison
#' @export
#'
#' @examples
#' # See vignette for complete example
compare_enrichment <- function(results_list,
                              padj_threshold = 0.05,
                              lfc_threshold = 1,
                              ...) {
  
  method_names <- names(results_list)
  enrichment_results <- list()
  
  # Run enrichment for each method
  for (method in method_names) {
    message(sprintf("\n=== Running enrichment for %s ===", method))
    
    enrichment_results[[method]] <- tryCatch({
      run_enrichment_from_de(
        results_list[[method]],
        padj_threshold = padj_threshold,
        lfc_threshold = lfc_threshold,
        ...
      )
    }, error = function(e) {
      warning(sprintf("Enrichment failed for %s: %s", method, e$message))
      NULL
    })
  }
  
  # Remove NULL entries
  enrichment_results <- enrichment_results[!sapply(enrichment_results, is.null)]
  
  # Find overlapping pathways
  if (length(enrichment_results) >= 2) {
    message("\n=== Computing pathway overlap ===")
    
    # Extract GO term IDs from each method
    go_terms <- lapply(enrichment_results, function(x) {
      if (!is.null(x$GO)) {
        x$GO$summary$ID
      } else {
        character(0)
      }
    })
    
    # Extract KEGG pathway IDs
    kegg_pathways <- lapply(enrichment_results, function(x) {
      if (!is.null(x$KEGG)) {
        x$KEGG$summary$ID
      } else {
        character(0)
      }
    })
    
    # Calculate overlap
    pathway_overlap <- list(
      GO = if (all(sapply(go_terms, length) > 0)) {
        compute_pathway_overlap(go_terms)
      } else NULL,
      KEGG = if (all(sapply(kegg_pathways, length) > 0)) {
        compute_pathway_overlap(kegg_pathways)
      } else NULL
    )
    
    enrichment_results$pathway_overlap <- pathway_overlap
  }
  
  return(enrichment_results)
}


#' Compute Pathway Overlap (Internal Helper)
#'
#' @param pathway_lists Named list of pathway ID vectors
#' @keywords internal
compute_pathway_overlap <- function(pathway_lists) {
  
  method_names <- names(pathway_lists)
  n_methods <- length(method_names)
  
  # Count occurrences
  all_pathways <- unique(unlist(pathway_lists))
  
  overlap_counts <- sapply(all_pathways, function(pw) {
    sum(sapply(pathway_lists, function(x) pw %in% x))
  })
  
  # Create summary
  consensus_pathways <- all_pathways[overlap_counts == n_methods]
  
  overlap_table <- table(overlap_counts)
  
  return(list(
    consensus_pathways = consensus_pathways,
    n_consensus = length(consensus_pathways),
    overlap_distribution = as.data.frame(overlap_table),
    all_pathways = all_pathways
  ))
}
