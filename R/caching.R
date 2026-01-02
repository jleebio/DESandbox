#' Cache Analysis Results
#'
#' Saves analysis results with metadata for reproducibility and faster reloading.
#'
#' @param results Analysis results (any R object)
#' @param cache_dir Directory to store cache files (default: "./desandbox_cache")
#' @param cache_name Optional name for cache file (default: auto-generated hash)
#' @param metadata Optional metadata list to store with results
#'
#' @return Path to cached file
#' @export
#'
#' @importFrom digest digest
#'
#' @examples
#' # See vignette for complete example
cache_analysis <- function(results,
                          cache_dir = "./desandbox_cache",
                          cache_name = NULL,
                          metadata = NULL) {
  
  # Create cache directory if it doesn't exist
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
    message(sprintf("Created cache directory: %s", cache_dir))
  }
  
  # Generate cache name if not provided
  if (is.null(cache_name)) {
    # Create hash from results
    results_hash <- digest::digest(results)
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    cache_name <- paste0("desandbox_", timestamp, "_", substr(results_hash, 1, 8))
  }
  
  # Add .rds extension if not present
  if (!grepl("\\.rds$", cache_name)) {
    cache_name <- paste0(cache_name, ".rds")
  }
  
  cache_path <- file.path(cache_dir, cache_name)
  
  # Prepare cache object
  cache_object <- list(
    results = results,
    cache_metadata = list(
      cached_date = Sys.time(),
      desandbox_version = tryCatch(
        as.character(utils::packageVersion("DESandbox")),
        error = function(e) "unknown"
      ),
      r_version = R.version.string,
      user_metadata = metadata
    )
  )
  
  # Save to file
  saveRDS(cache_object, cache_path)
  message(sprintf("Results cached to: %s", cache_path))
  
  return(invisible(cache_path))
}


#' Load Cached Analysis
#'
#' Loads previously cached analysis results.
#'
#' @param cache_path Path to cached .rds file
#' @param check_version Whether to check package version compatibility (default: TRUE)
#'
#' @return The cached results object
#' @export
#'
#' @examples
#' # See vignette for complete example
load_cached_analysis <- function(cache_path, check_version = TRUE) {
  
  if (!file.exists(cache_path)) {
    stop(sprintf("Cache file not found: %s", cache_path))
  }
  
  message(sprintf("Loading cache from: %s", cache_path))
  cache_object <- readRDS(cache_path)
  
  # Check structure
  if (!is.list(cache_object) || !"results" %in% names(cache_object)) {
    warning("Cache file may not be a valid DESandbox cache. Returning raw content.")
    return(cache_object)
  }
  
  # Display metadata
  if ("cache_metadata" %in% names(cache_object)) {
    meta <- cache_object$cache_metadata
    message(sprintf("Cached on: %s", meta$cached_date))
    message(sprintf("DESandbox version: %s", meta$desandbox_version))
    
    if (check_version) {
      current_version <- tryCatch(
        as.character(utils::packageVersion("DESandbox")),
        error = function(e) "unknown"
      )
      
      if (current_version != meta$desandbox_version && 
          meta$desandbox_version != "unknown") {
        warning(sprintf(
          "Cache was created with DESandbox v%s, current version is v%s. Results may be incompatible.",
          meta$desandbox_version, current_version
        ))
      }
    }
  }
  
  return(cache_object$results)
}


#' List Available Cache Files
#'
#' Lists all cache files in the specified directory.
#'
#' @param cache_dir Directory containing cache files (default: "./desandbox_cache")
#'
#' @return A data.frame with cache file information
#' @export
#'
#' @examples
#' list_cache_files()
list_cache_files <- function(cache_dir = "./desandbox_cache") {
  
  if (!dir.exists(cache_dir)) {
    message(sprintf("Cache directory does not exist: %s", cache_dir))
    return(data.frame())
  }
  
  cache_files <- list.files(cache_dir, pattern = "\\.rds$", full.names = TRUE)
  
  if (length(cache_files) == 0) {
    message("No cache files found")
    return(data.frame())
  }
  
  # Get file info
  file_info <- file.info(cache_files)
  
  cache_info <- data.frame(
    file = basename(cache_files),
    path = cache_files,
    size_mb = round(file_info$size / 1024^2, 2),
    modified = file_info$mtime,
    stringsAsFactors = FALSE
  )
  
  cache_info <- cache_info[order(cache_info$modified, decreasing = TRUE), ]
  
  return(cache_info)
}


#' Clear Cache
#'
#' Removes cache files older than specified age.
#'
#' @param cache_dir Directory containing cache files (default: "./desandbox_cache")
#' @param older_than Remove files older than this many days (default: NULL, removes all)
#' @param confirm Whether to ask for confirmation (default: TRUE)
#'
#' @return Number of files removed
#' @export
#'
#' @examples
#' # Clear cache files older than 30 days
#' clear_cache(older_than = 30)
clear_cache <- function(cache_dir = "./desandbox_cache",
                       older_than = NULL,
                       confirm = TRUE) {
  
  if (!dir.exists(cache_dir)) {
    message(sprintf("Cache directory does not exist: %s", cache_dir))
    return(0)
  }
  
  cache_files <- list.files(cache_dir, pattern = "\\.rds$", full.names = TRUE)
  
  if (length(cache_files) == 0) {
    message("No cache files to remove")
    return(0)
  }
  
  # Filter by age if specified
  if (!is.null(older_than)) {
    file_info <- file.info(cache_files)
    cutoff_date <- Sys.time() - (older_than * 24 * 60 * 60)
    cache_files <- cache_files[file_info$mtime < cutoff_date]
    
    if (length(cache_files) == 0) {
      message(sprintf("No cache files older than %d days", older_than))
      return(0)
    }
  }
  
  # Confirm deletion
  if (confirm) {
    response <- readline(prompt = sprintf(
      "Remove %d cache file(s)? (yes/no): ", length(cache_files)
    ))
    
    if (tolower(response) != "yes") {
      message("Operation cancelled")
      return(0)
    }
  }
  
  # Remove files
  n_removed <- sum(file.remove(cache_files))
  message(sprintf("Removed %d cache file(s)", n_removed))
  
  return(n_removed)
}


#' Generate Analysis Report Hash
#'
#' Creates a reproducible hash for tracking analysis parameters.
#'
#' @param desandbox_object The input data object
#' @param parameters List of analysis parameters
#'
#' @return A hash string
#' @export
#'
#' @importFrom digest digest
#'
#' @examples
#' # Internal use
generate_analysis_hash <- function(desandbox_object, parameters) {
  
  # Extract key information
  hash_components <- list(
    n_genes = nrow(desandbox_object),
    n_samples = ncol(desandbox_object),
    condition_col = S4Vectors::metadata(desandbox_object)$condition_column,
    parameters = parameters
  )
  
  analysis_hash <- digest::digest(hash_components)
  
  return(analysis_hash)
}
