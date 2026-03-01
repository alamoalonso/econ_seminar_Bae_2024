#!/usr/bin/env Rscript

#' Entry Point Script for Forecasting Workflow
#'
#' This script runs the complete forecasting workflow from start to finish.
#' It can be run from the command line or sourced from another R session.
#'
#' Usage:
#'   Rscript scripts/run_workflow.R [dataset]
#'
#' Arguments:
#'   dataset: "US" (default) or "EURO"
#'
#' Environment variables:
#'   DATASET: Override dataset selection (US or EURO)
#'   DEBUG: Set to "true" to enable debug logging
#'   SERIES_SUBSET: Comma-separated list of series to forecast (default: all)

# ============================================================================
# Setup
# ============================================================================

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

# Source all R files
source_dir <- function(path) {
  files <- list.files(path, pattern = "\\.R$", full.names = TRUE)
  for (f in files) {
    source(f)
  }
}

# Get project root (assuming script is in scripts/ subdirectory)
script_dir <- dirname(sys.frame(1)$ofile)
if (length(script_dir) == 0 || script_dir == "") {
  # Fallback if running interactively
  script_dir <- getwd()
  if (basename(script_dir) == "scripts") {
    script_dir <- dirname(script_dir)
  }
} else {
  script_dir <- dirname(script_dir)
}

setwd(script_dir)

# Source all R modules
cat("Loading R modules...\n")
source_dir("R")

# ============================================================================
# Parse arguments and configuration
# ============================================================================

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Determine dataset
dataset_arg <- if (length(args) > 0) args[1] else "US"
dataset_env <- Sys.getenv("DATASET", unset = "")
dataset_choice <- if (dataset_env != "") dataset_env else dataset_arg

# Create configuration
if (toupper(dataset_choice) == "EURO") {
  cat("Using Euro Area configuration\n")
  config <- config_euro_default()
} else {
  cat("Using US FRED-MD configuration\n")
  config <- config_us_default()
}

# Apply environment variable overrides
if (Sys.getenv("DEBUG", unset = "") == "true") {
  config$debug <- TRUE
  cat("Debug mode enabled\n")
}

# Parse series subset if provided
series_subset_env <- Sys.getenv("SERIES_SUBSET", unset = "")
if (series_subset_env != "") {
  config$series_list <- strsplit(series_subset_env, ",")[[1]]
  cat(sprintf("Forecasting subset of %d series\n", length(config$series_list)))
}

# ============================================================================
# Run workflow
# ============================================================================

cat("\n")
cat("========================================\n")
cat("Starting Forecasting Workflow\n")
cat("========================================\n")

# Run main workflow
results <- run_workflow(config)

# Compute unified evaluation (RMSE + tests) - OPTIMIZED
# This computes forecasts once and derives both RMSE and tests from same forecasts
cat("\n")
cat("Computing evaluation (RMSE + tests)...\n")
evaluation <- compute_evaluation(results, config)
rmse_results <- evaluation$rmse_results
tests_results <- evaluation$tests_results

# Compute MCS evaluation if enabled and forecasts are available
mcs_results <- NULL
if (isTRUE(config$mcs$enabled) && !is.null(evaluation$forecasts)) {
  cat("\n")
  cat("Computing MCS evaluation...\n")
  tryCatch({
    mcs_results <- compute_mcs_evaluation(evaluation$forecasts, config)
  }, error = function(e) {
    log_warn(sprintf("MCS evaluation failed: %s", e$message), config)
  })
}

# Save results
cat("\n")
cat("Saving results...\n")
output_dir <- save_results(results, rmse_results, config, tests_results,
                           forecasts = evaluation$forecasts,
                           mcs_results = mcs_results)

# Generate plots
cat("\n")
cat("Generating plots...\n")
tryCatch({
  generate_plots(rmse_results, config)
}, error = function(e) {
  log_warn(sprintf("Plot generation failed: %s", e$message), config)
})

# Print summary
print_summary(results, rmse_results, config, tests_results)

cat("\n")
cat(sprintf("Results saved to: %s\n", output_dir))
cat("\n")

# ============================================================================
# Exit
# ============================================================================

cat("Workflow complete!\n")
