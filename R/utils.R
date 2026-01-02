#' Global Variables and Imports
#'
#' This file defines global variables used in NSE (non-standard evaluation)
#' and imports commonly used functions from base R packages.
#'
#' @name utils
#' @keywords internal
NULL

# Import base R functions
#' @importFrom stats as.formula median setNames
#' @importFrom utils head
NULL

# Define global variables for NSE (used in ggplot2 and dplyr)
utils::globalVariables(c(
  # Visualization variables (ggplot2 aesthetics)
  "log2FC",
  "neg_log10_padj",
  "significance",
  "label",
  "x",
  "y",
  "set",
  "size",
  "Method",
  "Freq",
  "Direction"
))
