#!/usr/bin/env Rscript

#' Example Usage Script
#'
#' This script demonstrates various ways to use the refactored forecasting framework.

# ============================================================================
# Setup
# ============================================================================

# Set working directory to project root
# (adjust path as needed if running from elsewhere)
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

cat("Example Usage of Factor Forecasting Framework\n")
cat("==============================================\n\n")

# ============================================================================
# Example 1: Basic workflow with OPTIMIZED unified evaluation
# ============================================================================

cat("Example 1: Basic US workflow (OPTIMIZED)\n")
cat("-----------------------------------------\n")

# Create default configuration
config <- config_us_default()

# Optional: Customize for faster demo
# config$series_list <- c("INDPRO", "CPIAUCSL")
# config$horizons <- c(1, 12)
# config$do_tests <- TRUE  # Enable forecast comparison tests

# Run complete workflow
cat("Running workflow...\n")
results <- run_workflow(config)

# Compute RMSE + tests in one pass (OPTIMIZED - computes forecasts only once)
cat("Computing evaluation (RMSE + tests)...\n")
evaluation <- compute_evaluation(results, config)
rmse_results <- evaluation$rmse_results
tests_results <- evaluation$tests_results

# Display summary
print_summary(results, rmse_results, config, tests_results)

# Save results (including tests_results.csv, table5_dm.csv, table5_cw.csv if category file provided)
output_dir <- save_results(results, rmse_results, config, tests_results)

cat("\n")

# ============================================================================
# Example 1b: LEGACY approach (less efficient - for backward compatibility)
# ============================================================================

cat("Example 1b: Legacy approach (computes forecasts twice)\n")
cat("-------------------------------------------------------\n")

# NOTE: This approach is less efficient than compute_evaluation()
# Only use if you need RMSE or tests separately

cat("Computing RMSE only (legacy)...\n")
# rmse_results_legacy <- compute_rmse(results, config)

cat("Computing tests only (legacy)...\n")
# tests_results_legacy <- compute_tests(results, config)

cat("Skipped for efficiency - use compute_evaluation() instead!\n")
cat("\n")

# ============================================================================
# Example 2: Custom configuration for a subset of series
# ============================================================================

cat("Example 2: Subset of series with custom settings\n")
cat("-------------------------------------------------\n")

# Start with defaults and customize
config_subset <- config_us_default()
config_subset$series_list <- c("INDPRO", "CPIAUCSL", "S.P.500")
config_subset$horizons <- c(1, 6)
config_subset$factor_methods <- c("PCA")  # Only PCA, not PLS
config_subset$debug <- FALSE

cat(sprintf("Forecasting %d series at horizons %s\n",
            length(config_subset$series_list),
            paste(config_subset$horizons, collapse = ", ")))

# Run workflow
results_subset <- run_workflow(config_subset)
rmse_subset <- compute_rmse(results_subset, config_subset)

cat(sprintf("Generated %d RMSE rows\n", nrow(rmse_subset)))
cat("\n")

# ============================================================================
# Example 3: Running individual components
# ============================================================================

cat("Example 3: Running individual components\n")
cat("-----------------------------------------\n")

# Step-by-step execution for more control
config_manual <- config_us_default()

# Step 1: Load data
cat("Step 1: Loading data...\n")
data <- load_dataset(config_manual)
cat(sprintf("  Loaded %d observations x %d variables\n",
            nrow(data$raw_data), ncol(data$raw_data) - 1))

# Step 2: Preprocess
cat("Step 2: Preprocessing...\n")
dataset <- preprocess_dataset(data, config_manual)
cat(sprintf("  Balanced panel: %d predictors\n",
            length(dataset$balanced_predictors)))

# Step 3: Construct targets
cat("Step 3: Constructing targets...\n")
dataset$targets_list <- construct_targets(dataset$panel_final, config_manual$horizons)
cat(sprintf("  Targets for %d series\n", length(dataset$targets_list)))

# Step 4: Extract factors at a specific origin
cat("Step 4: Extracting factors at t=120...\n")
t_origin <- dataset$dates[120]
fcts <- extract_factors_at_origin(
  panel_final = dataset$panel_final,
  targets_list = dataset$targets_list,
  target_name = names(dataset$targets_list)[1],
  h = 1,
  t_origin = t_origin,
  k_max_pca = 12,
  k_max_pls = 12,
  config = config_manual
)
cat(sprintf("  PCA factors: %d x %d\n", nrow(fcts$pca$F), ncol(fcts$pca$F)))
cat(sprintf("  PLS factors: %d x %d\n", nrow(fcts$pls$F), ncol(fcts$pls$F)))

# Step 5: Run forecasts for a single series
cat("Step 5: Forecasting a single series...\n")
series_name <- names(dataset$targets_list)[1]
forecasts <- run_forecasts_for_series(
  series_name = series_name,
  h = 1,
  dates = dataset$dates,
  panel_final = dataset$panel_final,
  panel_std1 = dataset$panel_std1,
  targets_list = dataset$targets_list,
  factor_method = "PCA",
  k_max_pca = 12,
  k_max_pls = 12,
  config = config_manual
)
cat(sprintf("  Forecast results for %s: %d time points\n",
            series_name, length(forecasts$truth)))

cat("\n")

# ============================================================================
# Example 4: Debugging and diagnostics
# ============================================================================

cat("Example 4: Debugging mode\n")
cat("-------------------------\n")

# Enable debug output
config_debug <- config_us_default()
config_debug$debug <- TRUE
config_debug$trace_origins <- c(60, 120, 180)
config_debug$series_list <- c("INDPRO")
config_debug$horizons <- c(1)
config_debug$factor_methods <- c("PLS")  # Debug PLS specifically

cat("Running with debug mode enabled (traces at t=60, 120, 180)...\n")
cat("(Debug output will appear below)\n\n")

results_debug <- run_workflow(config_debug)
rmse_debug <- compute_rmse(results_debug, config_debug)

cat("\n")

# ============================================================================
# Example 5: Analyzing results
# ============================================================================

cat("Example 5: Analyzing RMSE results\n")
cat("----------------------------------\n")

# Use results from Example 1 (or recompute if needed)
if (exists("rmse_results")) {
  # Best performing model by horizon
  best_by_horizon <- rmse_results %>%
    filter(model != "AR") %>%  # Exclude benchmark
    group_by(h, scheme) %>%
    arrange(rmse_rel) %>%
    slice(1) %>%
    select(h, scheme, model, factor_method, k, rmse_rel)

  cat("Best performing models by horizon:\n")
  print(best_by_horizon, n = 20)

  # Mean RMSE across series for DIAR-LAG
  if (any(rmse_results$model == "DIAR-LAG")) {
    mean_diarlag <- rmse_results %>%
      filter(model == "DIAR-LAG", scheme == "recursive", h == 1) %>%
      group_by(factor_method, k) %>%
      summarise(mean_rmse_rel = mean(rmse_rel, na.rm = TRUE), .groups = "drop")

    cat("\nMean RMSE (relative to AR) for DIAR-LAG, h=1, recursive:\n")
    print(mean_diarlag)
  }
}

cat("\n")

# ============================================================================
# Example 6: Model comparison
# ============================================================================

cat("Example 6: Comparing PCA vs PLS\n")
cat("--------------------------------\n")

if (exists("rmse_results")) {
  pca_pls_comparison <- rmse_results %>%
    filter(model == "DIAR", k == 4, h == 1, scheme == "recursive") %>%
    group_by(factor_method) %>%
    summarise(
      mean_rmse = mean(rmse_rel, na.rm = TRUE),
      median_rmse = median(rmse_rel, na.rm = TRUE),
      n_series = n(),
      .groups = "drop"
    )

  cat("DIAR model (k=4, h=1, recursive) - PCA vs PLS:\n")
  print(pca_pls_comparison)
}

cat("\n")

# ============================================================================
# Summary
# ============================================================================

cat("Examples complete!\n")
cat("==================\n")
cat("You can modify these examples to:\n")
cat("- Test different series or horizons\n")
cat("- Compare different models\n")
cat("- Analyze performance patterns\n")
cat("- Debug specific issues\n")
cat("\n")
cat("See README.md for more detailed documentation.\n")
