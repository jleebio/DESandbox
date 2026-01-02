library(testthat)
library(DESandbox)

test_that("validate_counts works correctly", {
  # Create test data
  counts <- matrix(rpois(1000, 10), nrow = 100, ncol = 10)
  colnames(counts) <- paste0("sample", 1:10)
  rownames(counts) <- paste0("gene", 1:100)
  
  # Test basic validation
  result <- validate_counts(counts)
  expect_true(result$valid)
  expect_equal(result$n_samples, 10)
  
  # Test filtering
  counts_low <- counts
  counts_low[1:50, ] <- 1
  result <- validate_counts(counts_low, min_counts = 10)
  expect_lte(result$n_genes, 100)
  
  # Test error handling
  expect_error(validate_counts("not a matrix"))
  expect_error(validate_counts(matrix(-1, 10, 10), allow_negative = FALSE))
})

test_that("validate_metadata works correctly", {
  # Create test metadata
  metadata <- data.frame(
    sample = paste0("sample", 1:10),
    condition = factor(rep(c("A", "B"), each = 5)),
    batch = factor(rep(1:2, times = 5))
  )
  rownames(metadata) <- metadata$sample
  
  # Test basic validation
  result <- validate_metadata(metadata, condition_column = "condition")
  expect_true(result$valid)
  expect_equal(result$n_conditions, 2)
  
  # Test error handling
  expect_error(validate_metadata(metadata, condition_column = "missing"))
  
  metadata_single <- metadata
  metadata_single$condition <- factor("A")
  expect_error(validate_metadata(metadata_single, condition_column = "condition"))
})

test_that("create_desandbox_object works correctly", {
  # Create test data
  counts <- matrix(rpois(1000, 20), nrow = 100, ncol = 10)
  colnames(counts) <- paste0("sample", 1:10)
  rownames(counts) <- paste0("gene", 1:100)
  
  metadata <- data.frame(
    sample = colnames(counts),
    condition = factor(rep(c("control", "treatment"), each = 5))
  )
  rownames(metadata) <- metadata$sample
  
  # Create object
  dso <- create_desandbox_object(counts, metadata, "condition")
  
  expect_s4_class(dso, "SummarizedExperiment")
  expect_equal(nrow(dso), nrow(counts))
  expect_equal(ncol(dso), ncol(counts))
})

test_that("filter_results works correctly", {
  # Create test results
  results <- data.frame(
    gene = paste0("gene", 1:100),
    baseMean = runif(100, 10, 1000),
    log2FoldChange = rnorm(100, 0, 2),
    pvalue = runif(100, 0, 1),
    padj = runif(100, 0, 1)
  )
  
  # Test filtering
  filtered <- filter_results(
    results,
    padj_threshold = 0.05,
    lfc_threshold = 1
  )
  
  expect_true(all(filtered$padj < 0.05))
  expect_true(all(abs(filtered$log2FoldChange) >= 1))
  
  # Test direction filtering
  up <- filter_results(results, padj_threshold = 1, lfc_threshold = 0, direction = "up")
  expect_true(all(up$log2FoldChange > 0))
  
  down <- filter_results(results, padj_threshold = 1, lfc_threshold = 0, direction = "down")
  expect_true(all(down$log2FoldChange < 0))
})

test_that("standardize_results works correctly", {
  # Create mock DE result
  de_result <- list(
    method = "DESeq2",
    results = data.frame(
      gene = paste0("gene", 1:50),
      baseMean = runif(50, 10, 100),
      log2FoldChange = rnorm(50),
      pvalue = runif(50),
      padj = runif(50)
    )
  )
  
  standardized <- standardize_results(de_result)
  
  expect_true("method" %in% colnames(standardized))
  expect_equal(standardized$method[1], "DESeq2")
  expect_true(all(c("gene", "log2FoldChange", "pvalue", "padj") %in% colnames(standardized)))
})

test_that("compute_overlap works correctly", {
  # Create mock results from multiple methods
  results_list <- list(
    DESeq2 = list(
      results = data.frame(
        gene = paste0("gene", 1:50),
        log2FoldChange = rnorm(50, 2),
        pvalue = runif(50, 0, 0.01),
        padj = runif(50, 0, 0.01)
      )
    ),
    edgeR = list(
      results = data.frame(
        gene = paste0("gene", c(1:30, 51:70)),
        log2FoldChange = rnorm(50, 2),
        pvalue = runif(50, 0, 0.01),
        padj = runif(50, 0, 0.01)
      )
    )
  )
  
  overlap <- compute_overlap(
    results_list,
    padj_threshold = 0.05,
    lfc_threshold = 1
  )
  
  expect_true("pairwise_overlaps" %in% names(overlap))
  expect_true("overlap_summary" %in% names(overlap))
  expect_gt(overlap$overlap_summary$n_overlap[1], 0)
})

test_that("cache_analysis and load_cached_analysis work", {
  # Create temporary cache directory
  cache_dir <- tempdir()
  
  # Create test results
  test_results <- list(
    method = "test",
    data = data.frame(x = 1:10, y = 11:20)
  )
  
  # Cache results
  cache_path <- cache_analysis(
    test_results,
    cache_dir = cache_dir,
    cache_name = "test_cache"
  )
  
  expect_true(file.exists(cache_path))
  
  # Load cached results
  loaded_results <- load_cached_analysis(cache_path, check_version = FALSE)
  expect_equal(loaded_results$method, "test")
  expect_equal(nrow(loaded_results$data), 10)
  
  # Clean up
  unlink(cache_path)
})
