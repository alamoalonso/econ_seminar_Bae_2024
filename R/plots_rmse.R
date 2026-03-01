#' Plotting Functions for RMSE Results
#'
#' Functions to create publication-quality plots replicating Bae (2024) Figures 1 and 2.
#'
#' @name plots_rmse
NULL

#' Summarize mean RMSE by k
#'
#' Computes mean RMSE across series for each combination of scheme, factor_method, model, h, and k.
#' Enforces constraints on k (e.g., DIAR-LAG uses k <= 4).
#'
#' @param rmse_tbl Tibble with RMSE results from compute_rmse()
#' @param metric Character: "rmse_rel" (relative to AR) or "rmse" (absolute)
#' @param models Character vector: models to include (default: c("DI", "DIAR", "DIAR-LAG"))
#' @param horizons Integer vector: horizons to include (default: c(1, 6, 12, 24))
#' @param k_max_diarlag Integer: max k for DIAR-LAG (default: 4)
#' @return Tibble with columns: scheme, factor_method, model, h, k, mean_rmse, n_series
#' @export
#' @importFrom dplyr filter mutate group_by summarise ungroup
summarise_mean_rmse_by_k <- function(rmse_tbl,
                                     metric = c("rmse_rel", "rmse"),
                                     models = c("DI", "DIAR", "DIAR-LAG"),
                                     horizons = c(1, 6, 12, 24),
                                     k_max_diarlag = 4) {

  metric <- match.arg(metric)

  # Validate required columns
  required_cols <- c("series", "h", "scheme", "factor_method", "model", "k", "mse", "mse_ar", "rmse_rel")
  missing_cols <- setdiff(required_cols, names(rmse_tbl))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in rmse_tbl: ", paste(missing_cols, collapse = ", "))
  }

  # Compute rmse if needed
  if (metric == "rmse" && !("rmse" %in% names(rmse_tbl))) {
    rmse_tbl <- rmse_tbl %>%
      dplyr::mutate(rmse = sqrt(mse))
  }

  # Filter to requested models and horizons
  summary_tbl <- rmse_tbl %>%
    dplyr::filter(
      model %in% models,
      h %in% horizons,
      !is.na(k)  # Exclude AR rows which have k=NA
    )

  # Apply k constraint for DIAR-LAG
  summary_tbl <- summary_tbl %>%
    dplyr::filter(
      !(model == "DIAR-LAG" & k > k_max_diarlag)
    )

  # Select the appropriate metric column
  metric_col <- if (metric == "rmse_rel") "rmse_rel" else "rmse"

  # Group and summarize
  summary_tbl <- summary_tbl %>%
    dplyr::group_by(scheme, factor_method, model, h, k) %>%
    dplyr::summarise(
      mean_rmse = mean(.data[[metric_col]], na.rm = TRUE),
      n_series = sum(!is.na(.data[[metric_col]])),
      .groups = "drop"
    )

  summary_tbl
}

#' Plot mean RMSE by k for a single factor method
#'
#' Creates a 2x3 panel plot (2 rows for schemes, 3 columns for models) showing
#' mean RMSE as a function of k, with separate lines for each horizon.
#' Replicates the structure of Bae (2024) Figures 1 and 2.
#'
#' Each model column (DI, DIAR, DIAR-LAG) has independent y-axis scales.
#' This prevents DI's typically large k=1 values from compressing the
#' DIAR and DIAR-LAG panels, allowing clear visualization of trends across models.
#'
#' @param rmse_tbl Tibble with RMSE results from compute_rmse()
#' @param factor_method Character: "PCA" or "PLS"
#' @param metric Character: "rmse_rel" (relative to AR) or "rmse" (absolute)
#' @param schemes Character vector: c("recursive", "rolling")
#' @param models Character vector: c("DI", "DIAR", "DIAR-LAG")
#' @param horizons Integer vector: c(1, 6, 12, 24)
#' @param save_path Character: optional file path to save PNG (NULL = don't save)
#' @param width Numeric: plot width in inches (default: 10)
#' @param height Numeric: plot height in inches (default: 9)
#' @param dpi Numeric: resolution for saved plot (default: 300)
#' @return A patchwork/ggplot object (invisibly if saved)
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_point facet_wrap labs theme_minimal theme element_text element_blank scale_color_grey scale_linetype_manual scale_shape_manual ggsave unit guides guide_legend margin
#' @importFrom dplyr filter mutate
#' @importFrom patchwork wrap_plots plot_annotation plot_layout
plot_mean_rmse_kfactor <- function(rmse_tbl,
                                   factor_method = c("PCA", "PLS"),
                                   metric = c("rmse_rel", "rmse"),
                                   schemes = c("recursive", "rolling"),
                                   models = c("DI", "DIAR", "DIAR-LAG"),
                                   horizons = c(1, 6, 12, 24),
                                   save_path = NULL,
                                   series = NULL,
                                   width = 10,
                                   height = 9,
                                   dpi = 300) {

  factor_method <- match.arg(factor_method)
  metric <- match.arg(metric)

  # Prepare data
  plot_data <- summarise_mean_rmse_by_k(
    rmse_tbl = rmse_tbl,
    metric = metric,
    models = models,
    horizons = horizons
  )

  # Filter to requested factor method and schemes
  plot_data <- plot_data %>%
    dplyr::filter(
      factor_method == !!factor_method,
      scheme %in% schemes
    )

  if (nrow(plot_data) == 0) {
    warning(sprintf("No data available for factor_method=%s after filtering", factor_method))
    return(invisible(NULL))
  }

  # Ensure factor ordering
  plot_data <- plot_data %>%
    dplyr::mutate(
      scheme_label = factor(scheme,
                           levels = c("recursive", "rolling"),
                           labels = c("Recursive Estimation", "Rolling Estimation")),
      model = factor(model, levels = c("DI", "DIAR", "DIAR-LAG")),
      h = factor(h, levels = sort(unique(horizons)))
    )

  # Y-axis label based on metric
  y_label <- if (metric == "rmse_rel") {
    "Mean RMSE (relative to AR)"
  } else {
    "Mean RMSE"
  }

  # Create main title
  if (is.null(series)){
    main_title <- sprintf("Mean RMSE of k-%s", factor_method)
  } else {
    main_title <- sprintf("Mean RMSE of k-%s for series %s", factor_method, series)
  }
  # Create separate plots for each scheme to allow independent y-scales per model
  plot_list <- lapply(c("recursive", "rolling"), function(sch) {
    sch_label <- ifelse(sch == "recursive", "Recursive Estimation", "Rolling Estimation")

    data_subset <- plot_data %>% dplyr::filter(scheme == sch)

    if (nrow(data_subset) == 0) return(NULL)

    p <- ggplot2::ggplot(data_subset, ggplot2::aes(x = k, y = mean_rmse,
                                                     group = h,
                                                     color = h,
                                                     linetype = h,
                                                     shape = h)) +
      ggplot2::geom_line(linewidth = 0.6) +
      ggplot2::geom_point(size = 2) +
      ggplot2::facet_wrap(
        ~ model,
        nrow = 1,
        scales = "free_y",  # Independent y-scale per model!
        strip.position = "top"
      ) +
      ggplot2::scale_color_grey(
        name = "Horizon h",
        start = 0.2,
        end = 0.8
      ) +
      ggplot2::scale_linetype_manual(
        name = "Horizon h",
        values = c("solid", "dashed", "dotted", "dotdash")
      ) +
      ggplot2::scale_shape_manual(
        name = "Horizon h",
        values = c(16, 17, 15, 18)  # circle, triangle, square, diamond
      ) +
      ggplot2::labs(
        title = sch_label,
        x = "Number of factors k",
        y = y_label
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
        strip.text = ggplot2::element_text(face = "bold", size = 10),
        legend.position = "none",  # Will add shared legend at bottom
        panel.grid.minor = ggplot2::element_blank(),
        panel.spacing = ggplot2::unit(0.8, "lines")
      )

    p
  })

  # Remove NULL plots (if any scheme had no data)
  plot_list <- Filter(Negate(is.null), plot_list)

  if (length(plot_list) == 0) {
    warning("No plots generated - no data available")
    return(invisible(NULL))
  }

  # Re-create plots with legend enabled for the last one only
  # This allows patchwork to automatically collect and place the legend
  for (i in seq_along(plot_list)) {
    if (i == length(plot_list)) {
      # Last plot: show legend at bottom
      plot_list[[i]] <- plot_list[[i]] +
        ggplot2::theme(
          legend.position = "bottom",
          legend.box = "horizontal"
        ) +
        ggplot2::guides(
          color = ggplot2::guide_legend(nrow = 1, override.aes = list(linewidth = 1)),
          linetype = ggplot2::guide_legend(nrow = 1),
          shape = ggplot2::guide_legend(nrow = 1)
        )
    }
  }

  # Combine plots vertically with automatic legend collection
  final_plot <- patchwork::wrap_plots(
    plot_list,
    ncol = 1,
    guides = "collect"  # Collect legends from all plots
  ) +
    patchwork::plot_annotation(
      title = main_title,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14, margin = ggplot2::margin(b = 10))
      )
    ) +
    patchwork::plot_layout(
      heights = if (length(plot_list) == 2) c(1, 1) else c(1)
    )

  # Save if path provided
  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, final_plot, width = width, height = height, dpi = dpi)
    message(sprintf("Plot saved to: %s", save_path))
    return(invisible(final_plot))
  }

  final_plot
}

#' Plot Bae (2024) Figure 1: Mean RMSE of k-PCA
#'
#' Convenience wrapper to generate Figure 1 from Bae (2024) showing
#' PCA-based forecasting performance across schemes and models.
#'
#' @param rmse_tbl Tibble with RMSE results from compute_rmse()
#' @param metric Character: "rmse_rel" (relative to AR) or "rmse" (absolute)
#' @param save_dir Character: directory to save the plot
#' @param filename Character: optional custom filename (default: NULL uses "fig1_mean_rmse_k_pca.png")
#' @return A ggplot object (invisibly)
#' @export
plot_bae_fig1_kpca <- function(rmse_tbl,
                               metric = "rmse_rel",
                               save_dir,
                               filename = NULL,
                               config) {
  
  # Ensure config values exist
  if (is.null(config$schemes)) {
    sch <- c("recursive", "rolling")
  } else {
    sch <- config$schemes
  }
  
  if (is.null(config$models)) {
    mod <- c("DI", "DIAR", "DIAR-LAG")
  } else {
    mod <- config$models
  }
  
  if (is.null(config$horizons)) {
    hor <- c(1, 6, 12, 24)
  } else {
    hor <- config$horizons
  }
  
  # Construct filename
  if (is.null(filename)) {
    filename <- "fig1_mean_rmse_k_pca.png"
  }

  save_path <- file.path(save_dir, filename)

  # Ensure directory exists
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }

  # Generate plot
  plot_mean_rmse_kfactor(
    rmse_tbl = rmse_tbl,
    factor_method = "PCA",
    metric = metric,
    schemes = sch,
    models = mod,
    horizons = hor,
    save_path = save_path,
    series = NULL,
    width = 10,
    height = 9,
    dpi = 300
  )
  
  # Generate sub-plots for all series
  if(is.null(config$series_list)){
    series <- unique(rmse_tbl$series)
  } else {
    series <- config$series_list
  }
  
  for (s in series) {
    filename_series <- sub("\\.png$", paste0("_", s, ".png"), filename)
    save_path_series <- file.path(save_dir, filename_series)

    plot_mean_rmse_kfactor(
      rmse_tbl = rmse_tbl %>% filter(series == s),
      factor_method = "PCA",
      metric = metric,
      schemes = sch,
      models = mod,
      horizons = hor,
      save_path = save_path_series,
      series = s,
      width = 10,
      height = 9,
      dpi = 300
    )
  }
}

#' Plot Bae (2024) Figure 2: Mean RMSE of k-PLS
#'
#' Convenience wrapper to generate Figure 2 from Bae (2024) showing
#' PLS-based forecasting performance across schemes and models.
#'
#' @param rmse_tbl Tibble with RMSE results from compute_rmse()
#' @param metric Character: "rmse_rel" (relative to AR) or "rmse" (absolute)
#' @param save_dir Character: directory to save the plot
#' @param filename Character: optional custom filename (default: NULL uses "fig2_mean_rmse_k_pls.png")
#' @return A ggplot object (invisibly)
#' @export
plot_bae_fig2_kpls <- function(rmse_tbl,
                               metric = "rmse_rel",
                               save_dir,
                               filename = NULL,
                               config) {
  
  # Ensure config values exist
  if (is.null(config$schemes)) {
    sch <- c("recursive", "rolling")
  } else {
    sch <- config$schemes
  }
  
  if (is.null(config$models)) {
    mod <- c("DI", "DIAR", "DIAR-LAG")
  } else {
    mod <- config$models
  }
  
  if (is.null(config$horizons)) {
    hor <- c(1, 6, 12, 24)
  } else {
    hor <- config$horizons
  }

  # Construct filename
  if (is.null(filename)) {
    filename <- "fig2_mean_rmse_k_pls.png"
  }

  save_path <- file.path(save_dir, filename)

  # Ensure directory exists
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }

  # Generate plot
  plot_mean_rmse_kfactor(
    rmse_tbl = rmse_tbl,
    factor_method = "PLS",
    metric = metric,
    schemes = sch,
    models = mod,
    horizons = hor,
    save_path = save_path,
    series = NULL,
    width = 10,
    height = 9,
    dpi = 300
  )
  
  # Generate sub-plots for all series
  if(is.null(config$series_list)){
    series <- unique(rmse_tbl$series)
  } else {
    series <- config$series_list
  }
  
  for (s in series) {
    filename_series <- sub("\\.png$", paste0("_", s, ".png"), filename)
    save_path_series <- file.path(save_dir, filename_series)

    plot_mean_rmse_kfactor(
      rmse_tbl = rmse_tbl %>% filter(series == s),
      factor_method = "PLS",
      metric = metric,
      schemes = sch,
      models = mod,
      horizons = hor,
      save_path = save_path_series,
      series = s,
      width = 10,
      height = 9,
      dpi = 300
    )
  }
}

#' Generate plots for an already-completed run
#'
#' Convenience function to generate Bae-style plots for a run that has already
#' been completed and saved. Loads rmse_results.csv from the run directory and
#' generates plots.
#'
#' @param run_dir Character: path to the run output directory (e.g., "outputs/20231214_153045")
#' @param metric Character: "rmse_rel" (relative to AR) or "rmse" (absolute)
#' @param output_dir Character: optional output directory (default: NULL uses run_dir)
#' @return List with paths to generated plots (invisibly)
#' @export
#' @importFrom readr read_csv
#' @examples
#' \dontrun{
#' # Generate plots for a specific run
#' generate_plots_for_run("outputs/20231214_153045")
#'
#' # Generate plots for latest run
#' runs <- list.dirs("outputs", recursive = FALSE)
#' latest_run <- runs[length(runs)]
#' generate_plots_for_run(latest_run)
#' }
generate_plots_for_run <- function(run_dir,
                                   metric = "rmse_rel",
                                   output_dir = NULL) {

  # Validate run directory exists
  if (!dir.exists(run_dir)) {
    stop(sprintf("Run directory not found: %s", run_dir))
  }

  # Check for rmse_results.csv and config.rds
  rmse_file <- file.path(run_dir, "rmse_results.csv")
  if (!file.exists(rmse_file)) {
    stop(sprintf("rmse_results.csv not found in %s", run_dir))
  }
  
  config_file <- file.path(run_dir, "config.rds")
  if (!file.exists(config_file)) {
    stop(sprintf("config.rds not found in %s", run_dir))
  }

  # Load RMSE results and config
  message(sprintf("Loading RMSE results and config from: %s", rmse_file))
  rmse_results <- readr::read_csv(rmse_file, show_col_types = FALSE)
  config_run <- readRDS(config_file)

  # Determine output directory
  if (is.null(output_dir)) {
    output_dir <- run_dir
  }

  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Generate plots
  plots_generated <- list()

  # Generate PCA plot if PCA data exists
  if ("PCA" %in% rmse_results$factor_method) {
    message("Generating Figure 1 (k-PCA)...")
    tryCatch({
      plot_bae_fig1_kpca(
        rmse_results,
        metric = metric,
        save_dir = output_dir,
        config = config_run
      )
      fig1_path <- file.path(output_dir, "fig1_mean_rmse_k_pca.png")
      plots_generated$fig1 <- fig1_path
      message(sprintf("  -> Saved: %s", fig1_path))
    }, error = function(e) {
      warning(sprintf("Failed to generate PCA plot: %s", e$message))
    })
  } else {
    message("No PCA results found, skipping Figure 1")
  }

  # Generate PLS plot if PLS data exists
  if ("PLS" %in% rmse_results$factor_method) {
    message("Generating Figure 2 (k-PLS)...")
    tryCatch({
      plot_bae_fig2_kpls(
        rmse_results,
        metric = metric,
        save_dir = output_dir,
        config = config_run
      )
      fig2_path <- file.path(output_dir, "fig2_mean_rmse_k_pls.png")
      plots_generated$fig2 <- fig2_path
      message(sprintf("  -> Saved: %s", fig2_path))
    }, error = function(e) {
      warning(sprintf("Failed to generate PLS plot: %s", e$message))
    })
  } else {
    message("No PLS results found, skipping Figure 2")
  }

  if (length(plots_generated) == 0) {
    warning("No plots were generated. Check that rmse_results.csv contains valid data.")
  } else {
    message(sprintf("\nSuccessfully generated plot(s) in: %s",
                    output_dir))
  }

  invisible(plots_generated)
}

#' Generate plots for multiple runs
#'
#' Batch generate plots for multiple completed runs.
#'
#' @param run_dirs Character vector: paths to run output directories
#' @param metric Character: "rmse_rel" (relative to AR) or "rmse" (absolute)
#' @return List of lists with paths to generated plots for each run (invisibly)
#' @export
#' @examples
#' \dontrun{
#' # Generate plots for all runs in outputs directory
#' runs <- list.dirs("outputs", recursive = FALSE)
#' generate_plots_for_multiple_runs(runs)
#'
#' # Generate plots for specific runs
#' generate_plots_for_multiple_runs(c(
#'   "outputs/20231214_153045",
#'   "outputs/20231215_091530"
#' ))
#' }
generate_plots_for_multiple_runs <- function(run_dirs,
                                             metric = "rmse_rel") {

  message(sprintf("Generating plots for %d run(s)...\n", length(run_dirs)))

  results <- lapply(seq_along(run_dirs), function(i) {
    run_dir <- run_dirs[i]
    message(sprintf("=== Run %d/%d: %s ===", i, length(run_dirs), basename(run_dir)))

    tryCatch({
      generate_plots_for_run(run_dir, metric = metric)
    }, error = function(e) {
      warning(sprintf("Failed to process %s: %s", run_dir, e$message))
      NULL
    })
  })

  message("\n=== Batch processing complete ===")

  invisible(results)
}
