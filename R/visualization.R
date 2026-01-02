#' DESandbox Visualization Functions
#'
#' This module provides publication-quality visualizations for differential
#' expression analysis results including volcano plots and Venn diagrams.
#'
#' @name visualization
NULL

# ============================================================================
# Helper Functions
# ============================================================================

#' Classify Genes by Significance and Direction
#'
#' Internal helper to assign significance categories to genes.
#'
#' @param log2FC Numeric vector of log2 fold changes
#' @param padj Numeric vector of adjusted p-values
#' @param log2FC_threshold Log2 fold change threshold (default: 1)
#' @param padj_threshold Adjusted p-value threshold (default: 0.05)
#'
#' @return Factor with levels: "Upregulated", "Downregulated", "Not significant"
#' @keywords internal
classify_significance <- function(log2FC, padj, 
                                  log2FC_threshold = 1, 
                                  padj_threshold = 0.05) {
  
  # Handle NA values
  sig_cat <- rep("Not significant", length(log2FC))
  
  # Upregulated: log2FC > threshold AND padj < threshold
  sig_cat[!is.na(log2FC) & !is.na(padj) & 
          log2FC > log2FC_threshold & 
          padj < padj_threshold] <- "Upregulated"
  
  # Downregulated: log2FC < -threshold AND padj < threshold
  sig_cat[!is.na(log2FC) & !is.na(padj) & 
          log2FC < -log2FC_threshold & 
          padj < padj_threshold] <- "Downregulated"
  
  # Return as ordered factor for consistent plotting
  factor(sig_cat, 
         levels = c("Upregulated", "Downregulated", "Not significant"))
}


#' Identify Top N Genes by Significance
#'
#' Internal helper to select top significant genes for labeling.
#'
#' @param results_df Data frame with gene_id, log2FC, padj columns
#' @param n Number of top genes to select
#' @param by Criteria for ranking: "padj" or "abs_log2FC" (default: "padj")
#'
#' @return Character vector of gene IDs
#' @keywords internal
get_top_genes <- function(results_df, n = 10, by = "padj") {
  
  if (n <= 0 || nrow(results_df) == 0) {
    return(character(0))
  }
  
  # Remove rows with NA values
  clean_df <- results_df[!is.na(results_df$padj) & !is.na(results_df$log2FC), ]
  
  if (nrow(clean_df) == 0) {
    return(character(0))
  }
  
  # Rank genes
  if (by == "padj") {
    # Sort by padj (ascending), then by absolute log2FC (descending)
    clean_df <- clean_df[order(clean_df$padj, -abs(clean_df$log2FC)), ]
  } else if (by == "abs_log2FC") {
    # Sort by absolute log2FC (descending), then by padj (ascending)
    clean_df <- clean_df[order(-abs(clean_df$log2FC), clean_df$padj), ]
  }
  
  # Return top N gene IDs
  head(clean_df$gene_id, min(n, nrow(clean_df)))
}


#' Extract Significant Genes by Method and Direction
#'
#' Internal helper to filter genes for Venn diagram construction.
#'
#' @param results_df Data frame with method, gene_id, log2FC, padj columns
#' @param log2FC_threshold Log2 fold change threshold
#' @param padj_threshold Adjusted p-value threshold
#' @param direction Direction filter: "up", "down", or "both"
#'
#' @return Named list of character vectors (gene sets by method)
#' @keywords internal
extract_gene_sets <- function(results_df, 
                              log2FC_threshold = 1, 
                              padj_threshold = 0.05,
                              direction = "both") {
  
  # Validate direction
  direction <- match.arg(direction, c("up", "down", "both"))
  
  # Get unique methods
  methods <- unique(results_df$method)
  
  # Initialize result list
  gene_sets <- list()
  
  for (method_name in methods) {
    # Subset by method
    method_df <- results_df[results_df$method == method_name, ]
    
    # Apply filters based on direction
    if (direction == "up") {
      sig_genes <- method_df$gene_id[
        !is.na(method_df$log2FC) & 
        !is.na(method_df$padj) &
        method_df$log2FC > log2FC_threshold & 
        method_df$padj < padj_threshold
      ]
    } else if (direction == "down") {
      sig_genes <- method_df$gene_id[
        !is.na(method_df$log2FC) & 
        !is.na(method_df$padj) &
        method_df$log2FC < -log2FC_threshold & 
        method_df$padj < padj_threshold
      ]
    } else {  # both
      sig_genes <- method_df$gene_id[
        !is.na(method_df$log2FC) & 
        !is.na(method_df$padj) &
        abs(method_df$log2FC) > log2FC_threshold & 
        method_df$padj < padj_threshold
      ]
    }
    
    gene_sets[[method_name]] <- unique(sig_genes)
  }
  
  return(gene_sets)
}


#' Calculate Summary Statistics for Volcano Plot
#'
#' Internal helper to compute counts for plot annotations.
#'
#' @param results_df Data frame with classification column
#'
#' @return Data frame with counts by category
#' @keywords internal
compute_volcano_stats <- function(results_df) {
  
  if (!"significance" %in% colnames(results_df)) {
    stop("results_df must have 'significance' column")
  }
  
  # Count genes in each category
  stats <- as.data.frame(table(results_df$significance))
  colnames(stats) <- c("Category", "Count")
  
  return(stats)
}


# ============================================================================
# Main Visualization Functions
# ============================================================================

#' Create Volcano Plot for Differential Expression Results
#'
#' Generates publication-quality volcano plots showing log2 fold change
#' vs. statistical significance. Points are colored by significance category.
#'
#' @param results_df Data frame containing DE results with columns:
#'   \code{gene_id}, \code{log2FC}, \code{padj}, and optionally \code{method}
#' @param log2FC_threshold Log2 fold change threshold for significance (default: 1)
#' @param padj_threshold Adjusted p-value threshold for significance (default: 0.05)
#' @param facet_by_method If TRUE and method column exists, create faceted plot (default: TRUE)
#' @param label_top_n Number of top significant genes to label (default: 10, use 0 for no labels)
#' @param label_by Criteria for selecting top genes: "padj" or "abs_log2FC" (default: "padj")
#' @param point_size Size of points (default: 1.5)
#' @param point_alpha Transparency of points (default: 0.6)
#' @param colors Named vector of colors for categories (optional)
#' @param title Plot title (default: "Volcano Plot")
#' @param subtitle Plot subtitle (optional)
#'
#' @return A ggplot2 object
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline
#'   scale_color_manual theme_bw theme labs facet_wrap annotate
#' @importFrom ggrepel geom_text_repel
#'
#' @examples
#' \dontrun{
#' # Basic volcano plot
#' p <- plot_volcano(results_df)
#' print(p)
#'
#' # Customize thresholds and labeling
#' p <- plot_volcano(
#'   results_df,
#'   log2FC_threshold = 1.5,
#'   padj_threshold = 0.01,
#'   label_top_n = 20
#' )
#'
#' # Faceted by method
#' p <- plot_volcano(results_df, facet_by_method = TRUE)
#' }
plot_volcano <- function(results_df,
                        log2FC_threshold = 1,
                        padj_threshold = 0.05,
                        facet_by_method = TRUE,
                        label_top_n = 10,
                        label_by = "padj",
                        point_size = 1.5,
                        point_alpha = 0.6,
                        colors = NULL,
                        title = "Volcano Plot",
                        subtitle = NULL) {
  
  # Input validation
  required_cols <- c("gene_id", "log2FC", "padj")
  missing_cols <- setdiff(required_cols, colnames(results_df))
  
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", 
                paste(missing_cols, collapse = ", ")))
  }
  
  # Check for gene_id column (might be called "gene" in some results)
  if (!"gene_id" %in% colnames(results_df) && "gene" %in% colnames(results_df)) {
    results_df$gene_id <- results_df$gene
  }
  
  # Remove rows with NA in critical columns
  plot_df <- results_df[!is.na(results_df$log2FC) & !is.na(results_df$padj), ]
  
  if (nrow(plot_df) == 0) {
    stop("No valid data points after removing NAs")
  }
  
  # Add -log10(padj) for y-axis
  plot_df$neg_log10_padj <- -log10(plot_df$padj)
  
  # Cap infinite values (padj = 0 cases)
  max_finite <- max(plot_df$neg_log10_padj[is.finite(plot_df$neg_log10_padj)])
  plot_df$neg_log10_padj[is.infinite(plot_df$neg_log10_padj)] <- max_finite * 1.1
  
  # Classify significance
  plot_df$significance <- classify_significance(
    plot_df$log2FC,
    plot_df$padj,
    log2FC_threshold,
    padj_threshold
  )
  
  # Define default colors if not provided (Nature Genetics style)
  if (is.null(colors)) {
    colors <- c(
      "Upregulated" = "#D62728",      # Nature red
      "Downregulated" = "#1F77B4",    # Nature blue
      "Not significant" = "#D3D3D3"   # Light gray
    )
  }
  
  # Identify genes to label
  label_genes <- character(0)
  if (label_top_n > 0) {
    # Only label significant genes
    sig_df <- plot_df[plot_df$significance != "Not significant", ]
    if (nrow(sig_df) > 0) {
      label_genes <- get_top_genes(sig_df, label_top_n, by = label_by)
    }
  }
  
  # Add label column
  plot_df$label <- ""
  plot_df$label[plot_df$gene_id %in% label_genes] <- 
    plot_df$gene_id[plot_df$gene_id %in% label_genes]
  
  # Create base plot
  p <- ggplot2::ggplot(plot_df, 
                       ggplot2::aes(x = log2FC, 
                                   y = neg_log10_padj, 
                                   color = significance)) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::scale_color_manual(values = colors, name = "Classification") +
    
    # Add threshold lines with labels
    ggplot2::geom_hline(yintercept = -log10(padj_threshold), 
                       linetype = "dashed", 
                       color = "gray20", 
                       linewidth = 0.6) +
    ggplot2::geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), 
                       linetype = "dashed", 
                       color = "gray20", 
                       linewidth = 0.6) +
    
    # Add threshold annotations
    ggplot2::annotate("text", 
                     x = Inf, 
                     y = -log10(padj_threshold), 
                     label = sprintf("padj = %.3g", padj_threshold),
                     hjust = 1.1, 
                     vjust = -0.5, 
                     size = 3, 
                     color = "gray20",
                     fontface = "italic") +
    ggplot2::annotate("text", 
                     x = log2FC_threshold, 
                     y = Inf, 
                     label = sprintf("log2FC = %.2f", log2FC_threshold),
                     hjust = 1.1, 
                     vjust = 1.2, 
                     size = 3, 
                     color = "gray20",
                     angle = 90,
                     fontface = "italic") +
    ggplot2::annotate("text", 
                     x = -log2FC_threshold, 
                     y = Inf, 
                     label = sprintf("log2FC = %.2f", -log2FC_threshold),
                     hjust = -0.1, 
                     vjust = 1.2, 
                     size = 3, 
                     color = "gray20",
                     angle = 90,
                     fontface = "italic") +
    
    # Labels and theme
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = expression(log[2]~"fold change"),
      y = expression(-log[10]~"(adjusted"~italic(P)*"-value)")
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "gray95", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", linewidth = 0.8),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.5),
      legend.position = "top",
      legend.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.text = ggplot2::element_text(size = 9),
      legend.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.key = ggplot2::element_rect(fill = "white", color = NA),
      plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 10, color = "gray30", hjust = 0),
      axis.title = ggplot2::element_text(size = 10, face = "bold"),
      axis.text = ggplot2::element_text(size = 9, color = "black"),
      strip.background = ggplot2::element_rect(fill = "gray95", color = "black"),
      strip.text = ggplot2::element_text(size = 10, face = "bold")
    )
  
  # Add gene labels if requested
  if (label_top_n > 0 && length(label_genes) > 0) {
    label_df <- plot_df[plot_df$label != "", ]
    
    p <- p + ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(label = label),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "gray50",
      segment.size = 0.3,
      show.legend = FALSE
    )
  }
  
  # Add faceting if requested and method column exists
  if (facet_by_method && "method" %in% colnames(plot_df)) {
    p <- p + ggplot2::facet_wrap(~ method, ncol = 3)
  }
  
  # Add summary statistics as caption with threshold info
  stats <- compute_volcano_stats(plot_df)
  n_up <- stats$Count[stats$Category == "Upregulated"]
  n_down <- stats$Count[stats$Category == "Downregulated"]
  n_ns <- stats$Count[stats$Category == "Not significant"]
  
  caption <- sprintf(
    "Significance criteria: |log2FC| > %.2f, adjusted P < %.3g. Upregulated: %s | Downregulated: %s | Not significant: %s",
    log2FC_threshold,
    padj_threshold,
    ifelse(length(n_up) > 0, format(n_up, big.mark = ","), "0"),
    ifelse(length(n_down) > 0, format(n_down, big.mark = ","), "0"),
    ifelse(length(n_ns) > 0, format(n_ns, big.mark = ","), "0")
  )
  
  p <- p + ggplot2::labs(caption = caption) +
    ggplot2::theme(
      plot.caption = ggplot2::element_text(size = 8, hjust = 0, color = "gray30")
    )
  
  return(p)
}


#' Create Venn Diagram for DE Gene Overlap
#'
#' Generates Venn diagrams showing overlap of differentially expressed genes
#' across multiple methods. Supports 2-way and 3-way comparisons.
#'
#' @param results_df Data frame containing DE results with columns:
#'   \code{method}, \code{gene_id}, \code{log2FC}, \code{padj}
#' @param log2FC_threshold Log2 fold change threshold for significance (default: 1)
#' @param padj_threshold Adjusted p-value threshold for significance (default: 0.05)
#' @param direction Direction to filter: "up" (upregulated), "down" (downregulated),
#'   or "both" (default: "both")
#' @param methods Vector of method names to include (default: NULL, uses all methods)
#' @param colors Vector of colors for sets (optional)
#' @param title Plot title (default: "DE Gene Overlap")
#' @param return_type Return "plot" for visualization or "sets" for gene lists (default: "plot")
#'
#' @return Either a ggplot object (return_type = "plot") or a list of gene sets
#'   (return_type = "sets")
#' @export
#'
#' @importFrom ggplot2 ggplot aes theme_void coord_fixed
#' @importFrom ggforce geom_circle
#'
#' @examples
#' \dontrun{
#' # 3-way Venn diagram for all DE genes
#' p <- plot_venn_de(results_df)
#' print(p)
#'
#' # Venn diagram for upregulated genes only
#' p <- plot_venn_de(results_df, direction = "up")
#'
#' # Get gene sets instead of plot
#' gene_sets <- plot_venn_de(results_df, return_type = "sets")
#' }
plot_venn_de <- function(results_df,
                        log2FC_threshold = 1,
                        padj_threshold = 0.05,
                        direction = "both",
                        methods = NULL,
                        colors = NULL,
                        title = "DE Gene Overlap",
                        return_type = "plot") {
  
  # Input validation
  required_cols <- c("method", "gene_id", "log2FC", "padj")
  missing_cols <- setdiff(required_cols, colnames(results_df))
  
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", 
                paste(missing_cols, collapse = ", ")))
  }
  
  # Filter methods if specified
  if (!is.null(methods)) {
    results_df <- results_df[results_df$method %in% methods, ]
  }
  
  # Extract gene sets
  gene_sets <- extract_gene_sets(
    results_df,
    log2FC_threshold,
    padj_threshold,
    direction
  )
  
  if (length(gene_sets) == 0) {
    stop("No gene sets found. Check your filtering criteria.")
  }
  
  # Return gene sets if requested
  if (return_type == "sets") {
    return(gene_sets)
  }
  
  # Create Venn diagram based on number of sets
  n_sets <- length(gene_sets)
  
  if (n_sets < 2 || n_sets > 3) {
    stop(sprintf("Venn diagrams support 2-3 sets. Found %d sets.", n_sets))
  }
  
  # Define default colors
  if (is.null(colors)) {
    colors <- c("#D62728", "#1F77B4", "#2CA02C")[1:n_sets]  # Nature Genetics palette
  }
  
  # Calculate overlaps
  venn_data <- calculate_venn_overlaps(gene_sets)
  
  # Create plot
  if (n_sets == 2) {
    p <- create_venn_2way(venn_data, colors, title, direction, log2FC_threshold, padj_threshold)
  } else {
    p <- create_venn_3way(venn_data, colors, title, direction, log2FC_threshold, padj_threshold)
  }
  
  return(p)
}


#' Calculate Venn Diagram Overlaps
#'
#' Internal helper to compute all set intersections.
#'
#' @param gene_sets Named list of character vectors (gene sets)
#'
#' @return List containing set names, sizes, and overlaps
#' @keywords internal
calculate_venn_overlaps <- function(gene_sets) {
  
  set_names <- names(gene_sets)
  n_sets <- length(gene_sets)
  
  result <- list(
    names = set_names,
    n_sets = n_sets
  )
  
  if (n_sets == 2) {
    A <- gene_sets[[1]]
    B <- gene_sets[[2]]
    
    result$only_A <- length(setdiff(A, B))
    result$only_B <- length(setdiff(B, A))
    result$both <- length(intersect(A, B))
    
    result$genes_only_A <- setdiff(A, B)
    result$genes_only_B <- setdiff(B, A)
    result$genes_both <- intersect(A, B)
    
  } else if (n_sets == 3) {
    A <- gene_sets[[1]]
    B <- gene_sets[[2]]
    C <- gene_sets[[3]]
    
    # Calculate all intersections
    AB <- intersect(A, B)
    AC <- intersect(A, C)
    BC <- intersect(B, C)
    ABC <- intersect(AB, C)
    
    # Calculate unique regions
    result$only_A <- length(setdiff(setdiff(A, B), C))
    result$only_B <- length(setdiff(setdiff(B, A), C))
    result$only_C <- length(setdiff(setdiff(C, A), B))
    result$AB_only <- length(setdiff(AB, C))
    result$AC_only <- length(setdiff(AC, B))
    result$BC_only <- length(setdiff(BC, A))
    result$ABC <- length(ABC)
    
    # Store gene lists
    result$genes_only_A <- setdiff(setdiff(A, B), C)
    result$genes_only_B <- setdiff(setdiff(B, A), C)
    result$genes_only_C <- setdiff(setdiff(C, A), B)
    result$genes_AB_only <- setdiff(AB, C)
    result$genes_AC_only <- setdiff(AC, B)
    result$genes_BC_only <- setdiff(BC, A)
    result$genes_ABC <- ABC
  }
  
  return(result)
}


#' Create 2-Way Venn Diagram
#'
#' Internal helper to create 2-set Venn diagram using ggplot2.
#'
#' @param venn_data List containing overlap data
#' @param colors Vector of colors
#' @param title Plot title
#' @param direction Gene direction label
#' @param log2FC_threshold Log2FC threshold used
#' @param padj_threshold Adjusted p-value threshold used
#'
#' @return ggplot object
#' @keywords internal
create_venn_2way <- function(venn_data, colors, title, direction, 
                             log2FC_threshold = 1, padj_threshold = 0.05) {
  
  # Define circle positions and radius
  r <- 1.5
  d <- 1.2  # Distance between centers
  
  circles <- data.frame(
    x = c(-d/2, d/2),
    y = c(0, 0),
    r = c(r, r),
    set = venn_data$names,
    color = colors[1:2]
  )
  
  # Calculate label positions
  labels <- data.frame(
    x = c(-d/2 - r/2, d/2 + r/2, 0),
    y = c(0, 0, 0),
    label = c(
      sprintf("%s\n%d", venn_data$names[1], venn_data$only_A),
      sprintf("%s\n%d", venn_data$names[2], venn_data$only_B),
      sprintf("%d", venn_data$both)
    ),
    size = c(4, 4, 5)
  )
  
  # Create plot
  p <- ggplot2::ggplot() +
    ggforce::geom_circle(
      data = circles,
      ggplot2::aes(x0 = x, y0 = y, r = r, fill = set),
      alpha = 0.3,
      color = "black",
      linewidth = 0.8,
      show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(values = setNames(colors[1:2], venn_data$names)) +
    ggplot2::geom_text(
      data = labels,
      ggplot2::aes(x = x, y = y, label = label, size = size),
      fontface = "bold",
      show.legend = FALSE
    ) +
    ggplot2::scale_size_identity() +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("Significance: |log2FC| > %.2f, adjusted P < %.3g | Direction: %s | Total unique: %s",
                        log2FC_threshold, padj_threshold, direction,
                        format(venn_data$only_A + venn_data$only_B + venn_data$both, big.mark = ",")),
      caption = sprintf("Overlap: %s genes (%s%%)", 
                       format(venn_data$both, big.mark = ","),
                       round(100 * venn_data$both / (venn_data$only_A + venn_data$only_B + venn_data$both), 1))
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9, color = "gray30"),
      plot.caption = ggplot2::element_text(hjust = 0.5, size = 8, color = "gray30", face = "italic"),
      legend.position = "none",
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
  
  return(p)
}


#' Create 3-Way Venn Diagram
#'
#' Internal helper to create 3-set Venn diagram using ggplot2.
#'
#' @param venn_data List containing overlap data
#' @param colors Vector of colors
#' @param title Plot title
#' @param direction Gene direction label
#' @param log2FC_threshold Log2FC threshold used
#' @param padj_threshold Adjusted p-value threshold used
#'
#' @return ggplot object
#' @keywords internal
create_venn_3way <- function(venn_data, colors, title, direction,
                             log2FC_threshold = 1, padj_threshold = 0.05) {
  
  # Define circle positions for 3-way Venn
  r <- 1.5
  angle_offset <- pi / 2  # Start at top
  
  circles <- data.frame(
    x = c(
      r * cos(angle_offset),
      r * cos(angle_offset + 2*pi/3),
      r * cos(angle_offset + 4*pi/3)
    ),
    y = c(
      r * sin(angle_offset),
      r * sin(angle_offset + 2*pi/3),
      r * sin(angle_offset + 4*pi/3)
    ),
    r = rep(r, 3),
    set = venn_data$names,
    color = colors[1:3]
  )
  
  # Calculate label positions
  # Unique regions (outside overlaps)
  offset <- 1.8
  labels <- data.frame(
    x = c(
      # Unique to each set
      offset * cos(angle_offset),
      offset * cos(angle_offset + 2*pi/3),
      offset * cos(angle_offset + 4*pi/3),
      # Pairwise overlaps
      0.6 * cos(angle_offset + pi/3),
      0.6 * cos(angle_offset + pi),
      0.6 * cos(angle_offset + 5*pi/3),
      # Three-way overlap (center)
      0
    ),
    y = c(
      # Unique to each set
      offset * sin(angle_offset),
      offset * sin(angle_offset + 2*pi/3),
      offset * sin(angle_offset + 4*pi/3),
      # Pairwise overlaps
      0.6 * sin(angle_offset + pi/3),
      0.6 * sin(angle_offset + pi),
      0.6 * sin(angle_offset + 5*pi/3),
      # Three-way overlap (center)
      0
    ),
    label = c(
      sprintf("%s\n%d", venn_data$names[1], venn_data$only_A),
      sprintf("%s\n%d", venn_data$names[2], venn_data$only_B),
      sprintf("%s\n%d", venn_data$names[3], venn_data$only_C),
      sprintf("%d", venn_data$AB_only),
      sprintf("%d", venn_data$AC_only),
      sprintf("%d", venn_data$BC_only),
      sprintf("%d", venn_data$ABC)
    ),
    size = c(4, 4, 4, 4, 4, 4, 5),
    fontface = c("bold", "bold", "bold", "plain", "plain", "plain", "bold")
  )
  
  # Total unique genes
  total_genes <- venn_data$only_A + venn_data$only_B + venn_data$only_C +
                venn_data$AB_only + venn_data$AC_only + venn_data$BC_only +
                venn_data$ABC
  
  # Create plot
  p <- ggplot2::ggplot() +
    ggforce::geom_circle(
      data = circles,
      ggplot2::aes(x0 = x, y0 = y, r = r, fill = set),
      alpha = 0.3,
      color = "black",
      linewidth = 0.8,
      show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(values = setNames(colors[1:3], venn_data$names)) +
    ggplot2::geom_text(
      data = labels,
      ggplot2::aes(x = x, y = y, label = label, size = size),
      fontface = labels$fontface,
      show.legend = FALSE
    ) +
    ggplot2::scale_size_identity() +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("Significance: |log2FC| > %.2f, adjusted P < %.3g | Direction: %s | Total unique: %s",
                        log2FC_threshold, padj_threshold, direction,
                        format(total_genes, big.mark = ",")),
      caption = sprintf("Three-way consensus: %s genes (%s%% of total)", 
                       format(venn_data$ABC, big.mark = ","),
                       round(100 * venn_data$ABC / total_genes, 1))
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9, color = "gray30"),
      plot.caption = ggplot2::element_text(hjust = 0.5, size = 8, color = "gray30", face = "italic"),
      legend.position = "none",
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
  
  return(p)
}


#' Create Comparison Plot of DE Gene Counts
#'
#' Generates a bar plot comparing the number of significant genes across methods.
#'
#' @param results_df Data frame containing DE results with columns:
#'   \code{method}, \code{gene_id}, \code{log2FC}, \code{padj}
#' @param log2FC_threshold Log2 fold change threshold (default: 1)
#' @param padj_threshold Adjusted p-value threshold (default: 0.05)
#' @param split_by_direction If TRUE, split bars by up/down regulation (default: TRUE)
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#' \dontrun{
#' p <- plot_de_counts(results_df)
#' print(p)
#' }
plot_de_counts <- function(results_df,
                          log2FC_threshold = 1,
                          padj_threshold = 0.05,
                          split_by_direction = TRUE) {
  
  # Classify significance for each gene
  results_df$significance <- classify_significance(
    results_df$log2FC,
    results_df$padj,
    log2FC_threshold,
    padj_threshold
  )
  
  # Count by method and significance
  if (split_by_direction) {
    count_df <- as.data.frame(table(
      Method = results_df$method,
      Direction = results_df$significance
    ))
    count_df <- count_df[count_df$Direction != "Not significant", ]
    
    p <- ggplot2::ggplot(count_df, 
                        ggplot2::aes(x = Method, y = Freq, fill = Direction)) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::scale_fill_manual(
        values = c("Upregulated" = "#E64B35", "Downregulated" = "#4DBBD5")
      ) +
      ggplot2::labs(
        title = "Significant DE Genes by Method",
        y = "Number of Genes",
        x = "Method"
      )
  } else {
    # Count only significant genes (up + down)
    sig_df <- results_df[results_df$significance != "Not significant", ]
    count_df <- as.data.frame(table(Method = sig_df$method))
    
    p <- ggplot2::ggplot(count_df, ggplot2::aes(x = Method, y = Freq)) +
      ggplot2::geom_bar(stat = "identity", fill = "#00A087") +
      ggplot2::labs(
        title = "Total Significant DE Genes by Method",
        y = "Number of Genes",
        x = "Method"
      )
  }
  
  p <- p +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      panel.grid.major.x = ggplot2::element_blank()
    )
  
  return(p)
}
