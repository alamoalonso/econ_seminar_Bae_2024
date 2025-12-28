#!/usr/bin/env Rscript

#' Validation Script for Unified Evaluation Optimization
#'
#' This script validates that compute_evaluation() produces identical results
#' to separate compute_rmse() + compute_tests() calls, while being faster.
#'
#' Usage:
#'   Rscript scripts/validate_optimization.R

# ============================================================================
# Setup
# ============================================================================

# Set working directory to project root
if (basename(getwd()) == "scripts") {
  setwd("..")
}

# Load required libraries
suppressPackageStartupMessages({
  library(fbi)
  library(pls)
  library(tidyverse)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(purrr)
  library(lubridate)
})

# Source all R modules
source_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (f in source_files) source(f)

cat("\n")
cat("========================================\n")
cat("Validation: Unified Evaluation Function\n")
cat("========================================\n\n")

# ============================================================================
# Create test configuration (small subset for speed)
# ============================================================================

cat("Step 1: Creating test configuration\n")
cat("------------------------------------\n")

config <- config_us_default()
config$series_list <- c("INDPRO", "CPIAUCSL", "HOUST")  # Just 3 series
config$horizons <- c(1, 12)  # Just 2 horizons
config$factor_methods <- c("PCA", "PLS")  # Both methods
config$do_tests <- TRUE
config$test_types <- c("DM", "CW")
config$test_alpha <- c(0.10, 0.05, 0.01)

cat(sprintf("  Series: %d\n", length(config$series_list)))
cat(sprintf("  Horizons: %s\n", paste(config$horizons, collapse = ", ")))
cat(sprintf("  Factor methods: %s\n", paste(config$factor_methods, collapse = ", ")))
cat(sprintf("  Total combinations: %d\n\n",
            length(config$series_list) * length(config$horizons) * length(config$factor_methods)))

# ============================================================================
# Run workflow (shared data preparation)
# ============================================================================

cat("Step 2: Running workflow (data preparation)\n")
cat("--------------------------------------------\n")

results <- run_workflow(config)
cat(sprintf("  Dataset prepared with %d observations\n\n", nrow(results$dataset$panel_final)))

# ============================================================================
# Method 1: OPTIMIZED unified evaluation
# ============================================================================

cat("Step 3: Running OPTIMIZED unified evaluation\n")
cat("---------------------------------------------\n")

time_start_optimized <- Sys.time()
evaluation_optimized <- compute_evaluation(results, config)
time_end_optimized <- Sys.time()
time_optimized <- as.numeric(difftime(time_end_optimized, time_start_optimized, units = "secs"))

rmse_optimized <- evaluation_optimized$rmse_results
tests_optimized <- evaluation_optimized$tests_results

cat(sprintf("  RMSE results: %d rows\n", nrow(rmse_optimized)))
cat(sprintf("  Test results: %d rows\n", nrow(tests_optimized)))
cat(sprintf("  Time elapsed: %.2f seconds\n\n", time_optimized))

# ============================================================================
# Method 2: LEGACY separate calls
# ============================================================================

cat("Step 4: Running LEGACY separate calls\n")
cat("--------------------------------------\n")

time_start_legacy <- Sys.time()
rmse_legacy <- compute_rmse(results, config)
tests_legacy <- compute_tests(results, config)
time_end_legacy <- Sys.time()
time_legacy <- as.numeric(difftime(time_end_legacy, time_start_legacy, units = "secs"))

cat(sprintf("  RMSE results: %d rows\n", nrow(rmse_legacy)))
cat(sprintf("  Test results: %d rows\n", nrow(tests_legacy)))
cat(sprintf("  Time elapsed: %.2f seconds\n\n", time_legacy))

# ============================================================================
# Validation: Compare results
# ============================================================================

cat("Step 5: Validating results are identical\n")
cat("-----------------------------------------\n")

# Sort both datasets to ensure identical row order
rmse_optimized_sorted <- rmse_optimized %>%
  arrange(series, h, scheme, factor_method, model, k)

rmse_legacy_sorted <- rmse_legacy %>%
  arrange(series, h, scheme, factor_method, model, k)

tests_optimized_sorted <- tests_optimized %>%
  arrange(series_id, h, scheme, factor_method, model_class, k, test_type)

tests_legacy_sorted <- tests_legacy %>%
  arrange(series_id, h, scheme, factor_method, model_class, k, test_type)

# Compare RMSE results
cat("  Comparing RMSE results...\n")
rmse_comparison <- all.equal(rmse_optimized_sorted, rmse_legacy_sorted, tolerance = 1e-10)
if (isTRUE(rmse_comparison)) {
  cat("    ✓ RMSE results are IDENTICAL\n")
} else {
  cat("    ✗ RMSE results differ:\n")
  print(rmse_comparison)
}

# Compare test results
cat("  Comparing test results...\n")
tests_comparison <- all.equal(tests_optimized_sorted, tests_legacy_sorted, tolerance = 1e-10)
if (isTRUE(tests_comparison)) {
  cat("    ✓ Test results are IDENTICAL\n")
} else {
  cat("    ✗ Test results differ:\n")
  print(tests_comparison)
}

cat("\n")

# ============================================================================
# Performance comparison
# ============================================================================

cat("Step 6: Performance comparison\n")
cat("------------------------------\n")

speedup <- time_legacy / time_optimized
time_saved <- time_legacy - time_optimized
pct_saved <- 100 * (1 - time_optimized / time_legacy)

cat(sprintf("  Optimized time:  %.2f seconds\n", time_optimized))
cat(sprintf("  Legacy time:     %.2f seconds\n", time_legacy))
cat(sprintf("  Time saved:      %.2f seconds (%.1f%% faster)\n", time_saved, pct_saved))
cat(sprintf("  Speedup factor:  %.2fx\n\n", speedup))

# ============================================================================
# Summary
# ============================================================================

cat("========================================\n")
cat("Validation Summary\n")
cat("========================================\n")

if (isTRUE(rmse_comparison) && isTRUE(tests_comparison)) {
  cat("✓ VALIDATION PASSED\n")
  cat("  - Results are numerically identical\n")
  cat(sprintf("  - Performance improved by %.1f%%\n", pct_saved))
  cat("\n")
  cat("The optimized compute_evaluation() function is working correctly!\n")
  cat("It produces identical results while computing forecasts only once.\n")
} else {
  cat("✗ VALIDATION FAILED\n")
  cat("  - Results differ between methods\n")
  cat("  - Investigation required\n")
}

cat("========================================\n\n")

# ============================================================================
# Extrapolation to full dataset
# ============================================================================

cat("Extrapolation to full dataset:\n")
cat("-------------------------------\n")

# Assume full dataset has ~120 series, 4 horizons, 2 factor methods
full_combinations <- 120 * 4 * 2
test_combinations <- length(config$series_list) * length(config$horizons) * length(config$factor_methods)
scaling_factor <- full_combinations / test_combinations

estimated_time_optimized_full <- time_optimized * scaling_factor
estimated_time_legacy_full <- time_legacy * scaling_factor
estimated_time_saved_full <- estimated_time_legacy_full - estimated_time_optimized_full

cat(sprintf("  Test run: %d combinations\n", test_combinations))
cat(sprintf("  Full run: ~%d combinations (%.1fx more)\n", full_combinations, scaling_factor))
cat("\n")
cat(sprintf("  Estimated optimized time: %.1f minutes\n", estimated_time_optimized_full / 60))
cat(sprintf("  Estimated legacy time:    %.1f minutes\n", estimated_time_legacy_full / 60))
cat(sprintf("  Estimated time saved:     %.1f minutes\n\n", estimated_time_saved_full / 60))

cat("Note: Actual times may vary based on series complexity and system load.\n")
cat("\n")
