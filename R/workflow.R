#' Workflow Orchestration
#'
#' Top-level functions to run the complete forecasting workflow.
#'
#' @name workflow
NULL

#' Run complete forecasting workflow
#'
#' This is the main entry point for the entire forecasting pipeline.
#' It loads data, preprocesses, extracts factors, runs forecasts, and returns results.
#'
#' @param config Configuration list (from config_us_default() or config_euro_default())
#' @return A list with components:
#'   - config: The configuration used
#'   - dataset: The preprocessed dataset (panel_final, panel_std1, targets_list, etc.)
#'   - rmse_results: Tibble with RMSE results (if compute_rmse = TRUE)
#'   - run_timestamp: Character timestamp of the run
#' @export
#' @examples
#' config <- config_us_default()
#' results <- run_workflow(config)
run_workflow <- function(config) {
  log_info(paste(rep("=", 60), collapse = ""), config)
  log_info("Starting forecasting workflow", config)
  log_info(sprintf("Dataset: %s", config$dataset_id), config)
  log_info(sprintf("Run ID: %s", config$run_id), config)
  log_info(paste(rep("=", 60), collapse = ""), config)

  # Validate configuration
  validate_config(config)

  # Step 1: Load dataset
  log_info("Step 1: Loading dataset", config)
  data <- load_dataset(config)

  # Step 2: Preprocess dataset
  log_info("Step 2: Preprocessing dataset", config)
  dataset <- preprocess_dataset(data, config)

  # Step 3: Construct targets
  log_info("Step 3: Constructing h-step ahead targets", config)
  dataset$targets_list <- construct_targets(dataset$panel_final, config$horizons)

  log_info(sprintf("Targets constructed for %d series", length(dataset$targets_list)), config)

  # Step 4: Run forecasts (via compute_rmse, which calls run_forecasts_for_series internally)
  # This is done in compute_rmse() function

  # Return results
  results <- list(
    config = config,
    dataset = dataset,
    run_timestamp = Sys.time()
  )

  log_info("Workflow data preparation complete", config)

  results
}

#' Save results to disk
#'
#' Saves RMSE tables, test results, Table 5 summaries, forecasts, MCS results,
#' and optional plots to the output directory.
#'
#' @param results List returned from run_workflow()
#' @param rmse_results Tibble with RMSE results
#' @param config Configuration list
#' @param tests_results Tibble with test results (optional)
#' @param forecasts Tibble with OOS forecasts at origin level for MCS analysis (optional)
#' @param mcs_results Tibble with MCS evaluation results (optional)
#' @export
#' @importFrom readr write_csv
save_results <- function(results, rmse_results, config, tests_results = NULL,
                         forecasts = NULL, mcs_results = NULL) {
  log_info("Saving results to disk", config)

  # Create output directory if needed
  output_dir <- file.path(config$output_dir, config$run_id)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Save RMSE results
  rmse_file <- file.path(output_dir, "rmse_results.csv")
  readr::write_csv(rmse_results, rmse_file)
  log_info(sprintf("RMSE results saved to: %s", rmse_file), config)

  # Save test results if provided
  if (!is.null(tests_results)) {
    tests_file <- file.path(output_dir, "tests_results.csv")
    readr::write_csv(tests_results, tests_file)
    log_info(sprintf("Test results saved to: %s", tests_file), config)

    # Generate Table 5 summaries if category mapping is available
    if (!is.null(config$category_mapping_file) && file.exists(config$category_mapping_file)) {
      tryCatch({
        log_info("Generating Table 5 summaries", config)

        # Generate DM table
        if ("DM" %in% config$test_types) {
          table5_dm <- summarise_table5(tests_results, config, test_type = "DM")
          table5_dm_file <- file.path(output_dir, "table5_dm.csv")
          readr::write_csv(table5_dm, table5_dm_file)
          log_info(sprintf("Table 5 (DM) saved to: %s", table5_dm_file), config)
        }

        # Generate CW table
        if ("CW" %in% config$test_types) {
          table5_cw <- summarise_table5(tests_results, config, test_type = "CW")
          table5_cw_file <- file.path(output_dir, "table5_cw.csv")
          readr::write_csv(table5_cw, table5_cw_file)
          log_info(sprintf("Table 5 (CW) saved to: %s", table5_cw_file), config)
        }
      }, error = function(e) {
        log_warn(sprintf("Failed to generate Table 5 summaries: %s", e$message), config)
      })
    } else {
      log_info("Skipping Table 5 generation: category_mapping_file not provided", config)
    }
  }

  # Save forecasts if provided (for MCS analysis)
  if (!is.null(forecasts) && nrow(forecasts) > 0) {
    tryCatch({
      if (requireNamespace("arrow", quietly = TRUE)) {
        # Partitioned Parquet (efficient for large data, easy to query)
        forecasts_dir <- file.path(output_dir, "forecasts")
        arrow::write_dataset(
          forecasts,
          path = forecasts_dir,
          format = "parquet",
          partitioning = c("series_id", "h")
        )
        log_info(sprintf("Forecasts saved to: %s (partitioned parquet, %d rows)",
                         forecasts_dir, nrow(forecasts)), config)
      } else {
        # Fallback to CSV if arrow not available
        forecasts_file <- file.path(output_dir, "forecasts_long.csv")
        readr::write_csv(forecasts, forecasts_file)
        log_info(sprintf("Forecasts saved to: %s (%d rows)",
                         forecasts_file, nrow(forecasts)), config)
      }
    }, error = function(e) {
      log_warn(sprintf("Failed to save forecasts: %s", e$message), config)
    })
  }

  # Save MCS results if provided
  if (!is.null(mcs_results) && nrow(mcs_results) > 0) {
    tryCatch({
      if (requireNamespace("arrow", quietly = TRUE)) {
        mcs_file <- file.path(output_dir, "mcs_results.parquet")
        arrow::write_parquet(mcs_results, mcs_file)
        log_info(sprintf("MCS results saved to: %s (%d rows)", mcs_file, nrow(mcs_results)), config)
      } else {
        mcs_file <- file.path(output_dir, "mcs_results.csv")
        readr::write_csv(mcs_results, mcs_file)
        log_info(sprintf("MCS results saved to: %s (%d rows)", mcs_file, nrow(mcs_results)), config)
      }
    }, error = function(e) {
      log_warn(sprintf("Failed to save MCS results: %s", e$message), config)
    })
  }

  # Save configuration
  config_file <- file.path(output_dir, "config.rds")
  saveRDS(config, config_file)
  log_info(sprintf("Configuration saved to: %s", config_file), config)

  # Save summary statistics
  summary_file <- file.path(output_dir, "summary.txt")
  sink(summary_file)
  cat("Forecasting Workflow Summary\n")
  cat("============================\n\n")
  cat(sprintf("Dataset: %s\n", config$dataset_id))
  cat(sprintf("Run ID: %s\n", config$run_id))
  cat(sprintf("Run timestamp: %s\n", results$run_timestamp))
  cat(sprintf("Sample period: %s to %s\n", config$sample_start, config$sample_end))
  cat(sprintf("Evaluation period: %s to %s\n", config$eval_start, config$eval_end))
  cat(sprintf("Horizons: %s\n", paste(config$horizons, collapse = ", ")))
  cat(sprintf("Factor methods: %s\n", paste(config$factor_methods, collapse = ", ")))
  cat(sprintf("Number of series: %d\n", length(results$dataset$balanced_predictors)))
  cat(sprintf("Number of observations: %d\n", nrow(results$dataset$panel_final)))
  cat(sprintf("RMSE results: %d rows\n", nrow(rmse_results)))
  if (!is.null(forecasts)) {
    cat(sprintf("Forecasts saved: %d rows\n", nrow(forecasts)))
  }
  if (!is.null(mcs_results)) {
    cat(sprintf("MCS results: %d rows\n", nrow(mcs_results)))
  }
  cat("\n")
  sink()
  log_info(sprintf("Summary saved to: %s", summary_file), config)

  # Generate Bae-style plots if requested
  if (!is.null(config$make_plots) && config$make_plots) {
    log_info("Generating Bae-style RMSE plots", config)

    # Check if both PCA and PLS are in the results
    if ("PCA" %in% rmse_results$factor_method) {
      tryCatch({
        plot_bae_fig1_kpca(rmse_results, metric = "rmse_rel", save_dir = output_dir, config = config)
      }, error = function(e) {
        log_warn(sprintf("Failed to generate PCA plot: %s", e$message), config)
      })
    }

    if ("PLS" %in% rmse_results$factor_method) {
      tryCatch({
        plot_bae_fig2_kpls(rmse_results, metric = "rmse_rel", save_dir = output_dir, config = config)
      }, error = function(e) {
        log_warn(sprintf("Failed to generate PLS plot: %s", e$message), config)
      })
    }
  }

  log_info(sprintf("All results saved to: %s", output_dir), config)

  invisible(output_dir)
}

#' Generate summary plots from RMSE results
#'
#' Creates diagnostic plots for the forecasting results.
#'
#' @param rmse_results Tibble with RMSE results
#' @param config Configuration list
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs theme_minimal ggsave
#' @importFrom dplyr filter group_by summarise
generate_plots <- function(rmse_results, config) {
  log_info("Generating summary plots", config)

  output_dir <- file.path(config$output_dir, config$run_id)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Plot 1: Mean RMSE relative to AR across k for DIAR-LAG, h=1, recursive
  fig_data <- rmse_results %>%
    dplyr::filter(model == "DIAR-LAG", scheme == "recursive", h == 1) %>%
    dplyr::group_by(factor_method, k) %>%
    dplyr::summarise(mean_rmse_rel = mean(rmse_rel, na.rm = TRUE), .groups = "drop")

  if (nrow(fig_data) > 0) {
    p <- ggplot2::ggplot(fig_data,
                         ggplot2::aes(x = k, y = mean_rmse_rel, group = factor_method)) +
      ggplot2::geom_line(ggplot2::aes(linetype = factor_method)) +
      ggplot2::geom_point(ggplot2::aes(shape = factor_method)) +
      ggplot2::labs(
        x = "Number of factors k",
        y = "Mean RMSE relative to AR(BIC)",
        linetype = "Factors",
        shape    = "Factors",
        title    = "DIAR-LAG, recursive, h = 1: PCA vs PLS"
      ) +
      ggplot2::theme_minimal()

    plot_file <- file.path(output_dir, "diarlag_h1_recursive.png")
    ggplot2::ggsave(plot_file, p, width = 8, height = 6)
    log_info(sprintf("Plot saved to: %s", plot_file), config)
  } else {
    log_warn("No data available for DIAR-LAG plot", config)
  }

  # Additional plots can be added here

  invisible(NULL)
}

#' Print workflow summary to console
#'
#' @param results List returned from run_workflow()
#' @param rmse_results Tibble with RMSE results
#' @param config Configuration list
#' @param tests_results Tibble with test results (optional)
#' @export
print_summary <- function(results, rmse_results, config, tests_results = NULL) {
  cat("\n")
  cat("========================================\n")
  cat("Forecasting Workflow Complete\n")
  cat("========================================\n")
  cat(sprintf("Dataset: %s\n", config$dataset_id))
  cat(sprintf("Run ID: %s\n", config$run_id))
  cat(sprintf("Series: %d\n", length(results$dataset$balanced_predictors)))
  cat(sprintf("Observations: %d\n", nrow(results$dataset$panel_final)))
  cat(sprintf("Horizons: %s\n", paste(config$horizons, collapse = ", ")))
  cat(sprintf("Factor methods: %s\n", paste(config$factor_methods, collapse = ", ")))
  cat(sprintf("RMSE results: %d rows\n", nrow(rmse_results)))

  if (!is.null(tests_results)) {
    cat(sprintf("Test results: %d rows\n", nrow(tests_results)))
    cat(sprintf("Test types: %s\n", paste(unique(tests_results$test_type), collapse = ", ")))
  }

  cat("========================================\n")
  cat("\n")

  # Show a sample of RMSE results
  cat("Sample RMSE results (first 10 rows):\n")
  print(head(rmse_results, 10))
  cat("\n")

  # Show a sample of test results if available
  if (!is.null(tests_results)) {
    cat("Sample test results (first 10 rows):\n")
    print(head(tests_results, 10))
    cat("\n")
  }
}
