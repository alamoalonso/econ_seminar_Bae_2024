#' Evaluation Functions
#'
#' Functions for computing RMSE and evaluating forecast performance.
#'
#' @name evaluation
NULL

#' Compute MSE
#'
#' @param truth Numeric vector of true values
#' @param pred Numeric vector of predictions
#' @return MSE value, or NA if no valid observations
calc_mse <- function(truth, pred) {
  ok <- is.finite(truth) & is.finite(pred)
  if (!any(ok)) return(NA_real_)
  mean((truth[ok] - pred[ok])^2)
}

#' Compute RMSE for a single scheme (recursive or rolling)
#'
#' @param res List returned from run_forecasts_for_series()
#' @param dates Vector of dates
#' @param scheme Character: "recursive" or "rolling"
#' @param eval_start Date object: start of evaluation period
#' @param eval_end Date object: end of evaluation period
#' @param k_max_factor Integer: maximum number of factors
#' @param factor_method Character: "PCA" or "PLS"
#' @return Tibble with columns: scheme, factor_method, model, k, mse, mse_ar, rmse_rel
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows filter
compute_rmse_for_scheme <- function(res,
                                    dates,
                                    scheme = c("recursive", "rolling"),
                                    eval_start,
                                    eval_end,
                                    k_max_factor,
                                    factor_method) {
  scheme <- match.arg(scheme)

  # Select predictions based on scheme
  if (scheme == "recursive") {
    ar_pred   <- res$AR_rec
    DI_mat    <- res$DI_rec
    DIAR_mat  <- res$DIAR_rec
    DLAG_mat  <- res$DLAG_rec
  } else {
    ar_pred   <- res$AR_roll
    DI_mat    <- res$DI_roll
    DIAR_mat  <- res$DIAR_roll
    DLAG_mat  <- res$DLAG_roll
  }

  # Evaluation window
  idx_eval   <- which(dates >= eval_start & dates <= eval_end)
  truth_eval <- res$truth[idx_eval]

  # Compute AR benchmark MSE
  mse_ar <- calc_mse(truth_eval, ar_pred[idx_eval])
  if (!is.finite(mse_ar)) {
    warning("MSE(AR) is NA/Inf in this evaluation window.")
    return(NULL)
  }

  # Compute MSE for DI models (all k)
  di_list <- lapply(1:k_max_factor, function(k) {
    mse_k <- calc_mse(truth_eval, DI_mat[idx_eval, k])
    tibble::tibble(
      scheme        = scheme,
      factor_method = factor_method,
      model         = "DI",
      k             = k,
      mse           = mse_k,
      mse_ar        = mse_ar,
      rmse_rel      = mse_k / mse_ar
    )
  })

  # Compute MSE for DIAR models (all k)
  diar_list <- lapply(1:k_max_factor, function(k) {
    mse_k <- calc_mse(truth_eval, DIAR_mat[idx_eval, k])
    tibble::tibble(
      scheme        = scheme,
      factor_method = factor_method,
      model         = "DIAR",
      k             = k,
      mse           = mse_k,
      mse_ar        = mse_ar,
      rmse_rel      = mse_k / mse_ar
    )
  })

  # Compute MSE for DIAR-LAG models (k <= 4)
  dlag_list <- lapply(1:min(4, k_max_factor), function(k) {
    mse_k <- calc_mse(truth_eval, DLAG_mat[idx_eval, k])
    tibble::tibble(
      scheme        = scheme,
      factor_method = factor_method,
      model         = "DIAR-LAG",
      k             = k,
      mse           = mse_k,
      mse_ar        = mse_ar,
      rmse_rel      = mse_k / mse_ar
    )
  })

  # AR row
  ar_row <- tibble::tibble(
    scheme        = scheme,
    factor_method = factor_method,
    model         = "AR",
    k             = NA_integer_,
    mse           = mse_ar,
    mse_ar        = mse_ar,
    rmse_rel      = 1
  )

  # Combine all results
  dplyr::bind_rows(di_list, diar_list, dlag_list, ar_row) %>%
    dplyr::filter(is.finite(mse))
}

#' Compute RMSE for a single series and horizon
#'
#' @param series_name Character string
#' @param h Integer: forecast horizon
#' @param dates Vector of dates
#' @param panel_final data.frame with balanced predictors
#' @param panel_std1 data.frame with all predictors
#' @param targets_list Nested list of targets
#' @param factor_method Character: "PCA" or "PLS"
#' @param eval_start Date object
#' @param eval_end Date object
#' @param k_max_pca Integer
#' @param k_max_pls Integer
#' @param config Configuration list
#' @return Tibble with RMSE results
#' @export
#' @importFrom dplyr bind_rows mutate
compute_rmse_for_series_h <- function(series_name,
                                      h,
                                      dates,
                                      panel_final,
                                      panel_std1,
                                      targets_list,
                                      factor_method = c("PCA", "PLS"),
                                      eval_start = as.Date("1970-01-01"),
                                      eval_end   = as.Date("2019-12-01"),
                                      k_max_pca  = 12,
                                      k_max_pls  = 12,
                                      config = NULL) {

  factor_method <- match.arg(factor_method)
  k_max_factor  <- if (factor_method == "PCA") k_max_pca else k_max_pls

  # Run forecasts
  res <- run_forecasts_for_series(
    series_name   = series_name,
    h             = h,
    dates         = dates,
    panel_final   = panel_final,
    panel_std1    = panel_std1,
    targets_list  = targets_list,
    factor_method = factor_method,
    k_max_pca     = k_max_pca,
    k_max_pls     = k_max_pls,
    config        = config
  )

  # Compute RMSE for both schemes
  rec_tbl  <- compute_rmse_for_scheme(
    res, dates, "recursive",
    eval_start, eval_end,
    k_max_factor = k_max_factor,
    factor_method = factor_method
  )

  roll_tbl <- compute_rmse_for_scheme(
    res, dates, "rolling",
    eval_start, eval_end,
    k_max_factor = k_max_factor,
    factor_method = factor_method
  )

  # Combine and add series/horizon metadata
  dplyr::bind_rows(rec_tbl, roll_tbl) %>%
    dplyr::mutate(
      series = series_name,
      h      = h,
      .before = 1
    )
}

#' Compute RMSE across multiple series
#'
#' Main function for computing RMSE across all series in the results object.
#'
#' **NOTE**: For optimal performance when computing both RMSE and tests,
#' use compute_evaluation() instead, which computes forecasts only once.
#' This function is maintained for backward compatibility and standalone RMSE computation.
#'
#' @param results List of forecast results (from run_workflow)
#' @param config Configuration list
#' @return Tibble with RMSE results for all series/horizons/methods
#' @export
#' @importFrom purrr map_dfr
compute_rmse <- function(results, config) {
  log_info("Computing RMSE for all series and horizons", config)

  # Extract dataset components from results
  dataset <- results$dataset

  # Determine which series to evaluate
  if (is.null(config$series_list)) {
    series_vec <- dataset$balanced_predictors
  } else {
    series_vec <- config$series_list
  }

  # Calculate total iterations for progress tracking
  total_iterations <- length(series_vec) * length(config$horizons) * length(config$factor_methods)
  current_iteration <- 0

  # Compute RMSE for all combinations
  rmse_all <- purrr::map_dfr(series_vec, function(s) {
    purrr::map_dfr(config$horizons, function(h) {
      purrr::map_dfr(config$factor_methods, function(fm) {
        # Update progress counter
        current_iteration <<- current_iteration + 1
        cat(sprintf("\rComputing RMSE: [%d/%d] Series: %s, h=%d, Method: %s",
                    current_iteration, total_iterations, s, h, fm))
        flush.console()

        compute_rmse_for_series_h(
          series_name   = s,
          h             = h,
          dates         = dataset$dates,
          panel_final   = dataset$panel_final,
          panel_std1    = dataset$panel_std1,
          targets_list  = dataset$targets_list,
          factor_method = fm,
          eval_start    = config$eval_start,
          eval_end      = config$eval_end,
          k_max_pca     = config$k_max_pca,
          k_max_pls     = config$k_max_pls,
          config        = config
        )
      })
    })
  })

  # Print newline after progress is complete
  cat("\n")
  log_info(sprintf("RMSE computation complete: %d rows", nrow(rmse_all)), config)

  rmse_all
}

# ============================================================================
# Forecast Comparison Tests (Diebold-Mariano and Clark-West)
# ============================================================================

#' Compute Newey-West HAC variance of the mean
#'
#' @param x Numeric vector (loss differential series)
#' @param L Integer: truncation lag for HAC (L = max(h-1, 0))
#' @return Variance of the mean, or NA if invalid
nw_hac_var_mean <- function(x, L) {
  # Remove NAs
  x <- x[!is.na(x)]
  T <- length(x)

  if (T == 0) return(NA_real_)

  x_bar <- mean(x)

  # Compute gamma_0 (variance)
  gamma_0 <- mean((x - x_bar)^2)

  # If L = 0, no autocorrelation correction
  if (L == 0) {
    return(gamma_0 / T)
  }

  # Compute autocovariances gamma_ell for ell = 1, ..., L
  omega_hat <- gamma_0

  for (ell in 1:L) {
    # Bartlett weight
    w_ell <- 1 - ell / (L + 1)

    # Autocovariance at lag ell
    idx <- (ell + 1):T
    if (length(idx) == 0) next

    gamma_ell <- mean((x[idx] - x_bar) * (x[idx - ell] - x_bar))

    # Add to omega_hat
    omega_hat <- omega_hat + 2 * w_ell * gamma_ell
  }

  # Variance of the mean
  var_mean <- omega_hat / T

  # Return NA if non-positive (should not happen in theory, but numerical issues)
  if (!is.finite(var_mean) || var_mean <= 0) {
    return(NA_real_)
  }

  var_mean
}

#' Perform one-sided forecast comparison test (DM or CW)
#'
#' @param yt Numeric vector: realized values in evaluation window
#' @param f0 Numeric vector: benchmark forecast (AR)
#' @param f1 Numeric vector: competing forecast
#' @param h Integer: forecast horizon (for HAC lag L = max(h-1, 0))
#' @param test_type Character: "DM" or "CW"
#' @return List with stat, p_value_one_sided, n_oos
#' @export
forecast_test_one <- function(yt, f0, f1, h, test_type = c("DM", "CW")) {
  test_type <- match.arg(test_type)

  # Align data: remove NAs
  ok <- is.finite(yt) & is.finite(f0) & is.finite(f1)
  yt <- yt[ok]
  f0 <- f0[ok]
  f1 <- f1[ok]

  n_oos <- length(yt)

  # Minimum sample size check
  if (n_oos < 20) {
    return(list(
      stat = NA_real_,
      p_value_one_sided = NA_real_,
      n_oos = n_oos
    ))
  }

  # Compute forecast errors
  e0 <- yt - f0
  e1 <- yt - f1

  # Compute loss differential
  if (test_type == "DM") {
    # Diebold-Mariano: d_t = e0^2 - e1^2
    d <- e0^2 - e1^2
  } else {
    # Clark-West: d_t = e0^2 - (e1^2 - delta_f^2)
    delta_f <- f0 - f1
    d <- e0^2 - e1^2 + delta_f^2
  }

  # Compute mean differential
  d_bar <- mean(d)

  # Compute HAC variance with L = max(h - 1, 0)
  L <- max(h - 1, 0)
  var_d_bar <- nw_hac_var_mean(d, L)

  # Check for invalid variance
  if (is.na(var_d_bar) || var_d_bar <= 0) {
    return(list(
      stat = NA_real_,
      p_value_one_sided = NA_real_,
      n_oos = n_oos
    ))
  }

  # Compute test statistic
  stat <- d_bar / sqrt(var_d_bar)

  # One-sided p-value (H1: competing model is more accurate)
  # Under H0, stat ~ N(0,1) asymptotically
  # H1: E[d] > 0, so reject when stat is large
  p_value <- 1 - pnorm(stat)

  list(
    stat = stat,
    p_value_one_sided = p_value,
    n_oos = n_oos
  )
}

#' Compute forecast comparison tests for one series, horizon, scheme, and factor method
#'
#' @param res List returned from run_forecasts_for_series()
#' @param dates Vector of dates
#' @param scheme Character: "recursive" or "rolling"
#' @param eval_start Date object: start of evaluation period
#' @param eval_end Date object: end of evaluation period
#' @param k_max_factor Integer: maximum number of factors
#' @param factor_method Character: "PCA" or "PLS"
#' @param h Integer: forecast horizon
#' @param test_types Character vector: test types to compute (e.g., c("DM", "CW"))
#' @param test_alpha Numeric vector: significance levels (e.g., c(0.10, 0.05, 0.01))
#' @return Tibble with test results
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
compute_tests_for_scheme <- function(res,
                                     dates,
                                     scheme = c("recursive", "rolling"),
                                     eval_start,
                                     eval_end,
                                     k_max_factor,
                                     factor_method,
                                     h,
                                     test_types = c("DM", "CW"),
                                     test_alpha = c(0.10, 0.05, 0.01)) {
  scheme <- match.arg(scheme)

  # Select predictions based on scheme
  if (scheme == "recursive") {
    ar_pred   <- res$AR_rec
    DI_mat    <- res$DI_rec
    DIAR_mat  <- res$DIAR_rec
    DLAG_mat  <- res$DLAG_rec
  } else {
    ar_pred   <- res$AR_roll
    DI_mat    <- res$DI_roll
    DIAR_mat  <- res$DIAR_roll
    DLAG_mat  <- res$DLAG_roll
  }

  # Evaluation window
  idx_eval   <- which(dates >= eval_start & dates <= eval_end)
  truth_eval <- res$truth[idx_eval]
  ar_eval    <- ar_pred[idx_eval]

  # Initialize results list
  results <- list()

  # Test DI models (all k)
  for (k in 1:k_max_factor) {
    f1_eval <- DI_mat[idx_eval, k]

    for (test_type in test_types) {
      test_res <- forecast_test_one(truth_eval, ar_eval, f1_eval, h, test_type)

      # Create rejection indicators for each alpha
      reject_cols <- list()
      for (alpha in test_alpha) {
        col_name <- paste0("reject_", gsub("\\.", "", sprintf("%.2f", alpha)))
        reject_cols[[col_name]] <- !is.na(test_res$p_value_one_sided) &&
                                    test_res$p_value_one_sided < alpha
      }

      row <- tibble::tibble(
        scheme        = scheme,
        factor_method = factor_method,
        model_class   = "DI",
        k             = k,
        h             = h,
        test_type     = test_type,
        stat          = test_res$stat,
        p_value_one_sided = test_res$p_value_one_sided,
        n_oos         = test_res$n_oos
      )

      # Add rejection columns
      for (col_name in names(reject_cols)) {
        row[[col_name]] <- reject_cols[[col_name]]
      }

      results[[length(results) + 1]] <- row
    }
  }

  # Test DIAR models (all k)
  for (k in 1:k_max_factor) {
    f1_eval <- DIAR_mat[idx_eval, k]

    for (test_type in test_types) {
      test_res <- forecast_test_one(truth_eval, ar_eval, f1_eval, h, test_type)

      reject_cols <- list()
      for (alpha in test_alpha) {
        col_name <- paste0("reject_", gsub("\\.", "", sprintf("%.2f", alpha)))
        reject_cols[[col_name]] <- !is.na(test_res$p_value_one_sided) &&
                                    test_res$p_value_one_sided < alpha
      }

      row <- tibble::tibble(
        scheme        = scheme,
        factor_method = factor_method,
        model_class   = "DIAR",
        k             = k,
        h             = h,
        test_type     = test_type,
        stat          = test_res$stat,
        p_value_one_sided = test_res$p_value_one_sided,
        n_oos         = test_res$n_oos
      )

      for (col_name in names(reject_cols)) {
        row[[col_name]] <- reject_cols[[col_name]]
      }

      results[[length(results) + 1]] <- row
    }
  }

  # Test DIAR-LAG models (k <= 4)
  for (k in 1:min(4, k_max_factor)) {
    f1_eval <- DLAG_mat[idx_eval, k]

    for (test_type in test_types) {
      test_res <- forecast_test_one(truth_eval, ar_eval, f1_eval, h, test_type)

      reject_cols <- list()
      for (alpha in test_alpha) {
        col_name <- paste0("reject_", gsub("\\.", "", sprintf("%.2f", alpha)))
        reject_cols[[col_name]] <- !is.na(test_res$p_value_one_sided) &&
                                    test_res$p_value_one_sided < alpha
      }

      row <- tibble::tibble(
        scheme        = scheme,
        factor_method = factor_method,
        model_class   = "DIAR-LAG",
        k             = k,
        h             = h,
        test_type     = test_type,
        stat          = test_res$stat,
        p_value_one_sided = test_res$p_value_one_sided,
        n_oos         = test_res$n_oos
      )

      for (col_name in names(reject_cols)) {
        row[[col_name]] <- reject_cols[[col_name]]
      }

      results[[length(results) + 1]] <- row
    }
  }

  dplyr::bind_rows(results)
}

#' Compute forecast comparison tests across multiple series
#'
#' Main function for computing tests across all series in the results object.
#'
#' **NOTE**: For optimal performance when computing both RMSE and tests,
#' use compute_evaluation() instead, which computes forecasts only once.
#' This function is maintained for backward compatibility and standalone test computation.
#'
#' @param results List of forecast results (from run_workflow)
#' @param config Configuration list
#' @return Tibble with test results for all series/horizons/methods
#' @export
#' @importFrom purrr map_dfr
compute_tests <- function(results, config) {
  log_info("Computing forecast comparison tests for all series and horizons", config)

  # Extract dataset components from results
  dataset <- results$dataset

  # Determine which series to evaluate
  if (is.null(config$series_list)) {
    series_vec <- dataset$balanced_predictors
  } else {
    series_vec <- config$series_list
  }

  # Get test parameters from config
  test_types <- config$test_types
  test_alpha <- config$test_alpha

  # Calculate total iterations for progress tracking
  total_iterations <- length(series_vec) * length(config$horizons) * length(config$factor_methods)
  current_iteration <- 0

  # Compute tests for all combinations
  tests_all <- purrr::map_dfr(series_vec, function(s) {
    purrr::map_dfr(config$horizons, function(h) {
      purrr::map_dfr(config$factor_methods, function(fm) {
        # Update progress counter
        current_iteration <<- current_iteration + 1
        cat(sprintf("\rComputing tests: [%d/%d] Series: %s, h=%d, Method: %s",
                    current_iteration, total_iterations, s, h, fm))
        flush.console()

        # Determine k_max based on factor method
        k_max_factor <- if (fm == "PCA") config$k_max_pca else config$k_max_pls

        # Run forecasts for this series
        res <- run_forecasts_for_series(
          series_name   = s,
          h             = h,
          dates         = dataset$dates,
          panel_final   = dataset$panel_final,
          panel_std1    = dataset$panel_std1,
          targets_list  = dataset$targets_list,
          factor_method = fm,
          k_max_pca     = config$k_max_pca,
          k_max_pls     = config$k_max_pls,
          config        = config
        )

        # Compute tests for both schemes
        rec_tbl  <- compute_tests_for_scheme(
          res, dataset$dates, "recursive",
          config$eval_start, config$eval_end,
          k_max_factor = k_max_factor,
          factor_method = fm,
          h = h,
          test_types = test_types,
          test_alpha = test_alpha
        )

        roll_tbl <- compute_tests_for_scheme(
          res, dataset$dates, "rolling",
          config$eval_start, config$eval_end,
          k_max_factor = k_max_factor,
          factor_method = fm,
          h = h,
          test_types = test_types,
          test_alpha = test_alpha
        )

        # Combine and add series metadata
        dplyr::bind_rows(rec_tbl, roll_tbl) %>%
          dplyr::mutate(
            series_id = s,
            .before = 1
          )
      })
    })
  })

  # Print newline after progress is complete
  cat("\n")
  log_info(sprintf("Test computation complete: %d rows", nrow(tests_all)), config)

  tests_all
}

#' Compute unified evaluation (RMSE + tests) with optimized forecast computation
#'
#' This function computes forecasts once and derives both RMSE and test results
#' from the same forecast objects, eliminating redundant computation.
#'
#' @param results List of forecast results (from run_workflow)
#' @param config Configuration list
#' @return List with components:
#'   - rmse_results: Tibble with RMSE results
#'   - tests_results: Tibble with test results (NULL if do_tests = FALSE)
#' @export
#' @importFrom purrr map_dfr
#' @importFrom dplyr bind_rows mutate
compute_evaluation <- function(results, config) {
  log_info("Computing unified evaluation (RMSE + tests) for all series and horizons", config)

  # Extract dataset components from results
  dataset <- results$dataset

  # Determine which series to evaluate
  if (is.null(config$series_list)) {
    series_vec <- dataset$balanced_predictors
  } else {
    series_vec <- config$series_list
  }

  # Get test parameters from config
  test_types <- if (!is.null(config$test_types)) config$test_types else c("DM", "CW")
  test_alpha <- if (!is.null(config$test_alpha)) config$test_alpha else c(0.10, 0.05, 0.01)
  do_tests <- if (!is.null(config$do_tests)) config$do_tests else TRUE

  # Calculate total iterations for progress tracking
  total_iterations <- length(series_vec) * length(config$horizons) * length(config$factor_methods)
  current_iteration <- 0

  # Initialize result storage
  rmse_list <- list()
  tests_list <- list()

  # Single unified loop - compute forecasts once per combination
  for (s in series_vec) {
    for (h in config$horizons) {
      for (fm in config$factor_methods) {
        # Update progress counter
        current_iteration <- current_iteration + 1
        cat(sprintf("\rEvaluating: [%d/%d] Series: %s, h=%d, Method: %s",
                    current_iteration, total_iterations, s, h, fm))
        flush.console()

        # Determine k_max based on factor method
        k_max_factor <- if (fm == "PCA") config$k_max_pca else config$k_max_pls

        # ========================================
        # COMPUTE FORECASTS ONCE
        # ========================================
        res <- run_forecasts_for_series(
          series_name   = s,
          h             = h,
          dates         = dataset$dates,
          panel_final   = dataset$panel_final,
          panel_std1    = dataset$panel_std1,
          targets_list  = dataset$targets_list,
          factor_method = fm,
          k_max_pca     = config$k_max_pca,
          k_max_pls     = config$k_max_pls,
          config        = config
        )

        # ========================================
        # COMPUTE RMSE FROM FORECASTS
        # ========================================
        rmse_rec  <- compute_rmse_for_scheme(
          res, dataset$dates, "recursive",
          config$eval_start, config$eval_end,
          k_max_factor = k_max_factor,
          factor_method = fm
        )

        rmse_roll <- compute_rmse_for_scheme(
          res, dataset$dates, "rolling",
          config$eval_start, config$eval_end,
          k_max_factor = k_max_factor,
          factor_method = fm
        )

        # Combine and add series/horizon metadata
        if (!is.null(rmse_rec) || !is.null(rmse_roll)) {
          rmse_combined <- dplyr::bind_rows(rmse_rec, rmse_roll) %>%
            dplyr::mutate(
              series = s,
              h      = h,
              .before = 1
            )
          rmse_list[[length(rmse_list) + 1]] <- rmse_combined
        }

        # ========================================
        # COMPUTE TESTS FROM SAME FORECASTS
        # ========================================
        if (do_tests) {
          tests_rec  <- compute_tests_for_scheme(
            res, dataset$dates, "recursive",
            config$eval_start, config$eval_end,
            k_max_factor = k_max_factor,
            factor_method = fm,
            h = h,
            test_types = test_types,
            test_alpha = test_alpha
          )

          tests_roll <- compute_tests_for_scheme(
            res, dataset$dates, "rolling",
            config$eval_start, config$eval_end,
            k_max_factor = k_max_factor,
            factor_method = fm,
            h = h,
            test_types = test_types,
            test_alpha = test_alpha
          )

          # Combine and add series metadata
          tests_combined <- dplyr::bind_rows(tests_rec, tests_roll) %>%
            dplyr::mutate(
              series_id = s,
              .before = 1
            )
          tests_list[[length(tests_list) + 1]] <- tests_combined
        }
      }
    }
  }

  # Print newline after progress is complete
  cat("\n")

  # Combine all results
  rmse_results <- dplyr::bind_rows(rmse_list)
  tests_results <- if (length(tests_list) > 0) dplyr::bind_rows(tests_list) else NULL

  log_info(sprintf("Unified evaluation complete: RMSE=%d rows, Tests=%d rows",
                   nrow(rmse_results),
                   if (!is.null(tests_results)) nrow(tests_results) else 0), config)

  list(
    rmse_results = rmse_results,
    tests_results = tests_results
  )
}

#' Build Table 5 summary from test results
#'
#' Reproduces Bae (2024) Table 5 showing frequency and percentage of rejections
#' by category for a specific factor method, k, and test type.
#'
#' @param tests_tbl Tibble with test results from compute_tests()
#' @param config Configuration list
#' @param test_type Character: "DM" or "CW"
#' @return Tibble with Table 5 format (one per scheme and alpha level)
#' @export
#' @importFrom dplyr filter inner_join group_by summarise ungroup mutate
#' @importFrom tidyr pivot_longer
#' @importFrom readr read_csv
summarise_table5 <- function(tests_tbl, config, test_type = c("DM", "CW")) {
  test_type <- match.arg(test_type)

  # Load category mapping
  if (is.null(config$category_mapping_file)) {
    stop("category_mapping_file is NULL in config. Please provide a CSV with series,category columns.")
  }

  if (!file.exists(config$category_mapping_file)) {
    stop(sprintf("category_mapping_file does not exist: %s", config$category_mapping_file))
  }

  log_info(sprintf("Loading category mapping from: %s", config$category_mapping_file), config)
  category_map <- readr::read_csv(config$category_mapping_file, show_col_types = FALSE)

  if (!all(c("series", "category") %in% names(category_map))) {
    stop("category_mapping_file must contain 'series' and 'category' columns")
  }

  # Filter tests for Table 5 settings
  tbl5_data <- tests_tbl %>%
    dplyr::filter(
      factor_method == config$table5_factor_method,
      k == config$table5_k,
      test_type == test_type,
      model_class %in% c("DI", "DIAR", "DIAR-LAG")
    )

  if (nrow(tbl5_data) == 0) {
    stop(sprintf(
      "No test results found for factor_method=%s, k=%d, test_type=%s",
      config$table5_factor_method, config$table5_k, test_type
    ))
  }

  # Join with category mapping
  tbl5_data <- tbl5_data %>%
    dplyr::inner_join(category_map, by = c("series_id" = "series"))

  # Check for missing categories
  missing_series <- setdiff(unique(tbl5_data$series_id), category_map$series)
  if (length(missing_series) > 0) {
    stop(sprintf(
      "The following series are missing from category_mapping_file: %s",
      paste(missing_series, collapse = ", ")
    ))
  }

  # Build summary tables for each scheme and alpha
  results_list <- list()

  for (sch in unique(tbl5_data$scheme)) {
    for (alpha in config$test_alpha) {
      reject_col <- paste0("reject_", gsub("\\.", "", sprintf("%.2f", alpha)))

      if (!reject_col %in% names(tbl5_data)) {
        warning(sprintf("Rejection column %s not found in test results", reject_col))
        next
      }

      # Filter by scheme
      sch_data <- tbl5_data %>%
        dplyr::filter(scheme == sch)

      # Compute category-level summaries
      cat_summary <- sch_data %>%
        dplyr::group_by(category) %>%
        dplyr::summarise(
          frequency = sum(.data[[reject_col]], na.rm = TRUE),
          total = dplyr::n(),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          percentage = 100 * frequency / total,
          scheme = sch,
          alpha = alpha,
          test_type = test_type
        )

      # Compute overall summary
      overall <- sch_data %>%
        dplyr::summarise(
          category = "Overall",
          frequency = sum(.data[[reject_col]], na.rm = TRUE),
          total = dplyr::n(),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          percentage = 100 * frequency / total,
          scheme = sch,
          alpha = alpha,
          test_type = test_type
        )

      # Combine
      combined <- dplyr::bind_rows(overall, cat_summary)
      results_list[[length(results_list) + 1]] <- combined
    }
  }

  dplyr::bind_rows(results_list)
}
