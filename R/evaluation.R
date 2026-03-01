#' Evaluation Functions
#'
#' Functions for computing RMSE and evaluating forecast performance.
#'
#' @name evaluation
NULL

# Null-coalescing operator (define if not already available)
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}

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

# ============================================================================
# Forecast Persistence for MCS Analysis
# ============================================================================

#' Extract OOS forecasts from run_forecasts_for_series() result
#'
#' Converts the matrix-based output to a tidy long-format data frame
#' suitable for later loss computation and MCS analysis.
#'
#' @param res List from run_forecasts_for_series()
#' @param dates Vector of Date objects
#' @param series_id Character: series name
#' @param h Integer: forecast horizon
#' @param factor_method Character: "PCA" or "PLS"
#' @param eval_start Date: start of evaluation window
#' @param eval_end Date: end of evaluation window
#' @param k_max_factor Integer: maximum k used
#' @param run_id Character: run identifier
#' @return Tibble with forecast rows (may be NULL if no valid forecasts)
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows mutate case_when
extract_forecasts_from_res <- function(res, dates, series_id, h,
                                        factor_method, eval_start, eval_end,
                                        k_max_factor, run_id) {

  # Determine valid indices (where truth exists and is non-NA)
  valid_idx <- which(!is.na(res$truth))

  # Intersect with evaluation window
  eval_idx <- which(dates >= eval_start & dates <= eval_end)
  idx_to_extract <- intersect(valid_idx, eval_idx)

  if (length(idx_to_extract) == 0) return(NULL)

  # Pre-allocate list for efficiency

  rows_list <- vector("list")

  for (scheme in c("recursive", "rolling")) {
    ar_pred <- if (scheme == "recursive") res$AR_rec else res$AR_roll

    # AR rows: factor_method = NA (will be deduped later)
    for (t_idx in idx_to_extract) {
      rows_list[[length(rows_list) + 1]] <- tibble::tibble(
        run_id = run_id,
        series_id = series_id,
        h = h,
        scheme = scheme,
        factor_method = NA_character_,
        model_class = "AR",
        k = NA_integer_,
        origin_index = t_idx,
        origin_date = dates[t_idx],
        target_index = t_idx + h,
        target_date = dates[t_idx + h],
        y_true = res$truth[t_idx],
        y_hat = ar_pred[t_idx]
      )
    }

    # DI and DIAR (1:k_max_factor)
    for (model in c("DI", "DIAR")) {
      mat <- switch(paste0(model, "_", scheme),
        "DI_recursive" = res$DI_rec,
        "DI_rolling" = res$DI_roll,
        "DIAR_recursive" = res$DIAR_rec,
        "DIAR_rolling" = res$DIAR_roll
      )

      for (k in 1:k_max_factor) {
        for (t_idx in idx_to_extract) {
          rows_list[[length(rows_list) + 1]] <- tibble::tibble(
            run_id = run_id,
            series_id = series_id,
            h = h,
            scheme = scheme,
            factor_method = factor_method,
            model_class = model,
            k = k,
            origin_index = t_idx,
            origin_date = dates[t_idx],
            target_index = t_idx + h,
            target_date = dates[t_idx + h],
            y_true = res$truth[t_idx],
            y_hat = mat[t_idx, k]
          )
        }
      }
    }

    # DIAR-LAG (1:min(4, k_max_factor))
    mat <- if (scheme == "recursive") res$DLAG_rec else res$DLAG_roll
    k_max_dlag <- min(4, k_max_factor)

    for (k in 1:k_max_dlag) {
      for (t_idx in idx_to_extract) {
        rows_list[[length(rows_list) + 1]] <- tibble::tibble(
          run_id = run_id,
          series_id = series_id,
          h = h,
          scheme = scheme,
          factor_method = factor_method,
          model_class = "DIAR-LAG",
          k = k,
          origin_index = t_idx,
          origin_date = dates[t_idx],
          target_index = t_idx + h,
          target_date = dates[t_idx + h],
          y_true = res$truth[t_idx],
          y_hat = mat[t_idx, k]
        )
      }
    }
  }

  # Bind all rows at once (efficient)
  result <- dplyr::bind_rows(rows_list)

  # Add method_id
  result <- result %>%
    dplyr::mutate(
      method_id = dplyr::case_when(
        model_class == "AR" ~ paste(scheme, "AR", sep = "_"),
        TRUE ~ paste(scheme, factor_method, model_class, paste0("k", k), sep = "_")
      )
    )

  result
}


#' Extract OOS forecasts from run_forecasts_for_spec() result
#'
#' Converts the spec-based output to a tidy long-format data frame
#' suitable for MCS analysis. Handles both grid and dynamic specs.
#'
#' @param res List from run_forecasts_for_spec()
#' @param dates Vector of Date objects
#' @param series_id Character: series name
#' @param h Integer: forecast horizon
#' @param eval_start Date: start of evaluation window
#' @param eval_end Date: end of evaluation window
#' @param run_id Character: run identifier
#' @param first_forecast_idx Integer: first valid forecast index
#' @return Tibble with forecast rows (may be NULL if no valid forecasts)
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows mutate case_when
#' @export
extract_forecasts_from_spec_res <- function(res, dates, series_id, h,
                                             eval_start, eval_end, run_id,
                                             first_forecast_idx = 60) {

  spec_id <- res$spec_id
  factor_method <- res$factor_method
  k_mode <- res$k_mode
  k_rule <- res$k_rule
  k_max <- res$k_max

  # Determine k_selection_rule string for output
  k_selection_rule <- if (k_mode == "grid") "fixed" else k_rule

  # Determine valid indices (where truth exists and is non-NA)
  valid_idx <- which(!is.na(res$truth))

  # Intersect with evaluation window
  eval_idx <- which(dates >= eval_start & dates <= eval_end)
  idx_to_extract <- intersect(valid_idx, eval_idx)

  if (length(idx_to_extract) == 0) return(NULL)

  # Pre-allocate list for efficiency
  rows_list <- vector("list")

  for (scheme in c("recursive", "rolling")) {
    ar_pred <- if (scheme == "recursive") res$AR_rec else res$AR_roll

    # AR rows: factor_spec_id = NA (will be deduped later)
    for (t_idx in idx_to_extract) {
      # Compute training window
      training_start <- first_forecast_idx
      training_end <- t_idx

      rows_list[[length(rows_list) + 1]] <- tibble::tibble(
        run_id = run_id,
        series_id = series_id,
        h = h,
        scheme = scheme,
        factor_method = NA_character_,
        factor_spec_id = NA_character_,
        k_mode = NA_character_,
        k_selection_rule = "fixed",
        model_class = "AR",
        k = NA_integer_,
        origin_index = t_idx,
        origin_date = dates[t_idx],
        target_index = t_idx + h,
        target_date = if ((t_idx + h) <= length(dates)) dates[t_idx + h] else NA,
        training_window_start = training_start,
        training_window_end = training_end,
        y_true = res$truth[t_idx],
        y_hat = ar_pred[t_idx]
      )
    }

    if (k_mode == "grid") {
      # Grid mode: extract forecasts for k = 1..k_max
      k_max_factor <- k_max

      # DI and DIAR
      for (model in c("DI", "DIAR")) {
        mat <- switch(paste0(model, "_", scheme),
          "DI_recursive" = res$DI_rec,
          "DI_rolling" = res$DI_roll,
          "DIAR_recursive" = res$DIAR_rec,
          "DIAR_rolling" = res$DIAR_roll
        )

        for (k in 1:k_max_factor) {
          for (t_idx in idx_to_extract) {
            training_start <- first_forecast_idx
            training_end <- t_idx

            rows_list[[length(rows_list) + 1]] <- tibble::tibble(
              run_id = run_id,
              series_id = series_id,
              h = h,
              scheme = scheme,
              factor_method = factor_method,
              factor_spec_id = spec_id,
              k_mode = k_mode,
              k_selection_rule = k_selection_rule,
              model_class = model,
              k = k,
              origin_index = t_idx,
              origin_date = dates[t_idx],
              target_index = t_idx + h,
              target_date = if ((t_idx + h) <= length(dates)) dates[t_idx + h] else NA,
              training_window_start = training_start,
              training_window_end = training_end,
              y_true = res$truth[t_idx],
              y_hat = mat[t_idx, k]
            )
          }
        }
      }

      # DIAR-LAG
      mat <- if (scheme == "recursive") res$DLAG_rec else res$DLAG_roll
      k_max_dlag <- min(4, k_max_factor)

      for (k in 1:k_max_dlag) {
        for (t_idx in idx_to_extract) {
          training_start <- first_forecast_idx
          training_end <- t_idx

          rows_list[[length(rows_list) + 1]] <- tibble::tibble(
            run_id = run_id,
            series_id = series_id,
            h = h,
            scheme = scheme,
            factor_method = factor_method,
            factor_spec_id = spec_id,
            k_mode = k_mode,
            k_selection_rule = k_selection_rule,
            model_class = "DIAR-LAG",
            k = k,
            origin_index = t_idx,
            origin_date = dates[t_idx],
            target_index = t_idx + h,
            target_date = if ((t_idx + h) <= length(dates)) dates[t_idx + h] else NA,
            training_window_start = training_start,
            training_window_end = training_end,
            y_true = res$truth[t_idx],
            y_hat = mat[t_idx, k]
          )
        }
      }

    } else {
      # Dynamic mode: single forecast per origin with k = k_hat
      for (model in c("DI", "DIAR", "DIAR-LAG")) {
        pred <- switch(paste0(model, "_", scheme),
          "DI_recursive" = res$DI_rec,
          "DI_rolling" = res$DI_roll,
          "DIAR_recursive" = res$DIAR_rec,
          "DIAR_rolling" = res$DIAR_roll,
          "DIAR-LAG_recursive" = res$DLAG_rec,
          "DIAR-LAG_rolling" = res$DLAG_roll
        )

        for (t_idx in idx_to_extract) {
          # For DIAR-LAG, use k_hat_dlag (computed with k_max = 4)
          # For DI/DIAR, use k_hat (computed with full k_max)
          if (model == "DIAR-LAG") {
            k_used <- res$k_hat_dlag[t_idx]
          } else {
            k_used <- res$k_hat[t_idx]
          }

          # Skip if k is NA or invalid
          if (is.na(k_used)) next

          training_start <- first_forecast_idx
          training_end <- t_idx

          rows_list[[length(rows_list) + 1]] <- tibble::tibble(
            run_id = run_id,
            series_id = series_id,
            h = h,
            scheme = scheme,
            factor_method = factor_method,
            factor_spec_id = spec_id,
            k_mode = k_mode,
            k_selection_rule = k_selection_rule,
            model_class = model,
            k = k_used,  # Actual k used at this origin
            origin_index = t_idx,
            origin_date = dates[t_idx],
            target_index = t_idx + h,
            target_date = if ((t_idx + h) <= length(dates)) dates[t_idx + h] else NA,
            training_window_start = training_start,
            training_window_end = training_end,
            y_true = res$truth[t_idx],
            y_hat = pred[t_idx]
          )
        }
      }
    }
  }

  # Bind all rows at once (efficient)
  result <- dplyr::bind_rows(rows_list)

  if (nrow(result) == 0) return(NULL)

  # Add method_id: canonical identifier
  # Format: {scheme}_{factor_spec_id}_{model_class}[_k{n}] for grid
  # Format: {scheme}_{factor_spec_id}_{model_class} for dynamic
  result <- result %>%
    dplyr::mutate(
      method_id = dplyr::case_when(
        model_class == "AR" ~ paste(scheme, "AR", sep = "_"),
        k_mode == "grid" ~ paste(scheme, factor_spec_id, model_class, paste0("k", k), sep = "_"),
        TRUE ~ paste(scheme, factor_spec_id, model_class, sep = "_")
      )
    )

  result
}


#' Compute unified evaluation with factor_specs support
#'
#' This function computes forecasts once and derives both RMSE and test results
#' from the same forecast objects. Supports both grid and dynamic factor specs.
#'
#' @param results List of forecast results (from run_workflow)
#' @param config Configuration list
#' @return List with components:
#'   - rmse_results: Tibble with RMSE results
#'   - tests_results: Tibble with test results (NULL if do_tests = FALSE)
#'   - forecasts: Tibble with OOS forecasts at origin level (NULL if save_forecasts = FALSE)
#' @export
#' @importFrom purrr map_dfr
#' @importFrom dplyr bind_rows mutate filter distinct
compute_evaluation_with_specs <- function(results, config) {
  log_info("Computing unified evaluation (RMSE + tests) with factor_specs support", config)

  # Extract dataset components from results
  dataset <- results$dataset

  # Get factor_specs (with backward compatibility)
  factor_specs <- get_factor_specs(config)

  # Determine which series to evaluate
  if (is.null(config$series_list)) {
    series_vec <- dataset$balanced_predictors
  } else {
    series_vec <- config$series_list
  }

  # Validate that all series exist in targets_list
  available_series <- names(dataset$targets_list)
  missing_series <- setdiff(series_vec, available_series)

  if (length(missing_series) > 0) {
    log_warn(sprintf(
      "Skipping %d missing series: %s",
      length(missing_series), paste(head(missing_series, 5), collapse = ", ")
    ), config)
    series_vec <- intersect(series_vec, available_series)
  }

  if (length(series_vec) == 0) {
    stop("No valid series to evaluate after filtering.")
  }

  # Get evaluation parameters
  test_types <- config$test_types %||% c("DM", "CW")
  test_alpha <- config$test_alpha %||% c(0.10, 0.05, 0.01)
  do_tests <- config$do_tests %||% TRUE
  do_save_forecasts <- isTRUE(config$save_forecasts)

  # Progress tracking
  total_iterations <- length(series_vec) * length(config$horizons) * length(factor_specs)
  current_iteration <- 0

  # Initialize storage
  rmse_list <- list()
  tests_list <- list()
  forecasts_list <- if (do_save_forecasts) list() else NULL

  for (s in series_vec) {
    for (h in config$horizons) {
      for (spec in factor_specs) {
        current_iteration <- current_iteration + 1
        cat(sprintf("\rEvaluating: [%d/%d] Series: %s, h=%d, Spec: %s",
                    current_iteration, total_iterations, s, h, spec$id))
        flush.console()

        # Run forecasts for this spec
        res <- run_forecasts_for_spec(
          series_name = s,
          h = h,
          dates = dataset$dates,
          panel_final = dataset$panel_final,
          panel_std1 = dataset$panel_std1,
          targets_list = dataset$targets_list,
          factor_spec = spec,
          config = config
        )

        # Compute RMSE for this spec
        rmse_chunk <- compute_rmse_for_spec(
          res = res,
          dates = dataset$dates,
          eval_start = config$eval_start,
          eval_end = config$eval_end,
          series_name = s,
          h = h
        )
        if (!is.null(rmse_chunk)) {
          rmse_list[[length(rmse_list) + 1]] <- rmse_chunk
        }

        # Compute tests if requested
        if (do_tests) {
          tests_chunk <- compute_tests_for_spec(
            res = res,
            dates = dataset$dates,
            eval_start = config$eval_start,
            eval_end = config$eval_end,
            series_name = s,
            h = h,
            test_types = test_types,
            test_alpha = test_alpha
          )
          if (!is.null(tests_chunk)) {
            tests_list[[length(tests_list) + 1]] <- tests_chunk
          }
        }

        # Extract forecasts for MCS
        if (do_save_forecasts) {
          forecast_chunk <- extract_forecasts_from_spec_res(
            res = res,
            dates = dataset$dates,
            series_id = s,
            h = h,
            eval_start = config$eval_start,
            eval_end = config$eval_end,
            run_id = config$run_id,
            first_forecast_idx = config$first_forecast_idx %||% 60
          )
          if (!is.null(forecast_chunk)) {
            forecasts_list[[length(forecasts_list) + 1]] <- forecast_chunk
          }
        }
      }
    }
  }

  cat("\n")

  # Combine results
  rmse_results <- dplyr::bind_rows(rmse_list)
  tests_results <- if (length(tests_list) > 0) dplyr::bind_rows(tests_list) else NULL

  # Process forecasts: bind and deduplicate AR rows
  forecasts_all <- NULL
  if (do_save_forecasts && length(forecasts_list) > 0) {
    forecasts_all <- dplyr::bind_rows(forecasts_list)

    # Deduplicate AR rows (AR is computed identically for each factor_spec)
    ar_rows <- forecasts_all %>%
      dplyr::filter(model_class == "AR")
    non_ar_rows <- forecasts_all %>%
      dplyr::filter(model_class != "AR")

    # Keep only unique AR rows
    ar_deduped <- ar_rows %>%
      dplyr::distinct(run_id, series_id, h, scheme, origin_index, .keep_all = TRUE)

    forecasts_all <- dplyr::bind_rows(ar_deduped, non_ar_rows)

    log_info(sprintf("Forecasts collected: %d rows (after AR deduplication)",
                     nrow(forecasts_all)), config)
  }

  log_info(sprintf("Evaluation complete: RMSE=%d rows, Tests=%d rows",
                   nrow(rmse_results),
                   if (!is.null(tests_results)) nrow(tests_results) else 0), config)

  list(
    rmse_results = rmse_results,
    tests_results = tests_results,
    forecasts = forecasts_all
  )
}


#' Compute RMSE for a spec-based result
#'
#' @param res List from run_forecasts_for_spec()
#' @param dates Vector of dates
#' @param eval_start Date: start of evaluation
#' @param eval_end Date: end of evaluation
#' @param series_name Character: series name
#' @param h Integer: horizon
#' @return Tibble with RMSE results
#' @keywords internal
compute_rmse_for_spec <- function(res, dates, eval_start, eval_end, series_name, h) {

  spec_id <- res$spec_id
  factor_method <- res$factor_method
  k_mode <- res$k_mode
  k_rule <- res$k_rule
  k_max <- res$k_max

  # Evaluation window
  idx_eval <- which(dates >= eval_start & dates <= eval_end)
  truth_eval <- res$truth[idx_eval]

  results <- list()

  for (scheme in c("recursive", "rolling")) {
    # AR benchmark
    ar_pred <- if (scheme == "recursive") res$AR_rec else res$AR_roll
    mse_ar <- calc_mse(truth_eval, ar_pred[idx_eval])

    if (!is.finite(mse_ar)) next

    # AR row
    results[[length(results) + 1]] <- tibble::tibble(
      series = series_name,
      h = h,
      scheme = scheme,
      factor_method = NA_character_,
      factor_spec_id = NA_character_,
      k_mode = NA_character_,
      k_selection_rule = "fixed",
      model = "AR",
      k = NA_integer_,
      mse = mse_ar,
      mse_ar = mse_ar,
      rmse_rel = 1
    )

    if (k_mode == "grid") {
      # Grid mode: compute for each k
      for (model in c("DI", "DIAR", "DIAR-LAG")) {
        k_max_model <- if (model == "DIAR-LAG") min(4, k_max) else k_max

        for (k in 1:k_max_model) {
          mat <- switch(paste0(model, "_", scheme),
            "DI_recursive" = res$DI_rec,
            "DI_rolling" = res$DI_roll,
            "DIAR_recursive" = res$DIAR_rec,
            "DIAR_rolling" = res$DIAR_roll,
            "DIAR-LAG_recursive" = res$DLAG_rec,
            "DIAR-LAG_rolling" = res$DLAG_roll
          )

          mse_k <- calc_mse(truth_eval, mat[idx_eval, k])

          results[[length(results) + 1]] <- tibble::tibble(
            series = series_name,
            h = h,
            scheme = scheme,
            factor_method = factor_method,
            factor_spec_id = spec_id,
            k_mode = k_mode,
            k_selection_rule = "fixed",
            model = model,
            k = k,
            mse = mse_k,
            mse_ar = mse_ar,
            rmse_rel = mse_k / mse_ar
          )
        }
      }
    } else {
      # Dynamic mode: aggregate MSE across all k_hat values
      for (model in c("DI", "DIAR", "DIAR-LAG")) {
        pred <- switch(paste0(model, "_", scheme),
          "DI_recursive" = res$DI_rec,
          "DI_rolling" = res$DI_roll,
          "DIAR_recursive" = res$DIAR_rec,
          "DIAR_rolling" = res$DIAR_roll,
          "DIAR-LAG_recursive" = res$DLAG_rec,
          "DIAR-LAG_rolling" = res$DLAG_roll
        )

        mse_dyn <- calc_mse(truth_eval, pred[idx_eval])

        results[[length(results) + 1]] <- tibble::tibble(
          series = series_name,
          h = h,
          scheme = scheme,
          factor_method = factor_method,
          factor_spec_id = spec_id,
          k_mode = k_mode,
          k_selection_rule = k_rule,
          model = model,
          k = NA_integer_,  # k varies by origin for dynamic
          mse = mse_dyn,
          mse_ar = mse_ar,
          rmse_rel = mse_dyn / mse_ar
        )
      }
    }
  }

  if (length(results) == 0) return(NULL)
  dplyr::bind_rows(results) %>% dplyr::filter(is.finite(mse))
}


#' Compute forecast tests for a spec-based result
#'
#' @param res List from run_forecasts_for_spec()
#' @param dates Vector of dates
#' @param eval_start Date: start of evaluation
#' @param eval_end Date: end of evaluation
#' @param series_name Character: series name
#' @param h Integer: horizon
#' @param test_types Character vector: test types to compute
#' @param test_alpha Numeric vector: significance levels
#' @return Tibble with test results
#' @keywords internal
compute_tests_for_spec <- function(res, dates, eval_start, eval_end,
                                    series_name, h, test_types, test_alpha) {

  spec_id <- res$spec_id
  factor_method <- res$factor_method
  k_mode <- res$k_mode
  k_rule <- res$k_rule
  k_max <- res$k_max

  # Evaluation window
  idx_eval <- which(dates >= eval_start & dates <= eval_end)
  truth_eval <- res$truth[idx_eval]

  results <- list()

  for (scheme in c("recursive", "rolling")) {
    ar_pred <- if (scheme == "recursive") res$AR_rec else res$AR_roll
    ar_eval <- ar_pred[idx_eval]

    if (k_mode == "grid") {
      # Grid mode: test each k
      for (model in c("DI", "DIAR", "DIAR-LAG")) {
        k_max_model <- if (model == "DIAR-LAG") min(4, k_max) else k_max

        for (k in 1:k_max_model) {
          mat <- switch(paste0(model, "_", scheme),
            "DI_recursive" = res$DI_rec,
            "DI_rolling" = res$DI_roll,
            "DIAR_recursive" = res$DIAR_rec,
            "DIAR_rolling" = res$DIAR_roll,
            "DIAR-LAG_recursive" = res$DLAG_rec,
            "DIAR-LAG_rolling" = res$DLAG_roll
          )

          f1_eval <- mat[idx_eval, k]

          for (test_type in test_types) {
            test_res <- forecast_test_one(truth_eval, ar_eval, f1_eval, h, test_type)

            reject_cols <- list()
            for (alpha in test_alpha) {
              col_name <- paste0("reject_", gsub("\\.", "", sprintf("%.2f", alpha)))
              reject_cols[[col_name]] <- !is.na(test_res$p_value_one_sided) &&
                                          test_res$p_value_one_sided < alpha
            }

            row <- tibble::tibble(
              series_id = series_name,
              scheme = scheme,
              factor_method = factor_method,
              factor_spec_id = spec_id,
              k_mode = k_mode,
              k_selection_rule = "fixed",
              model_class = model,
              k = k,
              h = h,
              test_type = test_type,
              stat = test_res$stat,
              p_value_one_sided = test_res$p_value_one_sided,
              n_oos = test_res$n_oos
            )

            for (col_name in names(reject_cols)) {
              row[[col_name]] <- reject_cols[[col_name]]
            }

            results[[length(results) + 1]] <- row
          }
        }
      }
    } else {
      # Dynamic mode: single test per model
      for (model in c("DI", "DIAR", "DIAR-LAG")) {
        pred <- switch(paste0(model, "_", scheme),
          "DI_recursive" = res$DI_rec,
          "DI_rolling" = res$DI_roll,
          "DIAR_recursive" = res$DIAR_rec,
          "DIAR_rolling" = res$DIAR_roll,
          "DIAR-LAG_recursive" = res$DLAG_rec,
          "DIAR-LAG_rolling" = res$DLAG_roll
        )

        f1_eval <- pred[idx_eval]

        for (test_type in test_types) {
          test_res <- forecast_test_one(truth_eval, ar_eval, f1_eval, h, test_type)

          reject_cols <- list()
          for (alpha in test_alpha) {
            col_name <- paste0("reject_", gsub("\\.", "", sprintf("%.2f", alpha)))
            reject_cols[[col_name]] <- !is.na(test_res$p_value_one_sided) &&
                                        test_res$p_value_one_sided < alpha
          }

          row <- tibble::tibble(
            series_id = series_name,
            scheme = scheme,
            factor_method = factor_method,
            factor_spec_id = spec_id,
            k_mode = k_mode,
            k_selection_rule = k_rule,
            model_class = model,
            k = NA_integer_,
            h = h,
            test_type = test_type,
            stat = test_res$stat,
            p_value_one_sided = test_res$p_value_one_sided,
            n_oos = test_res$n_oos
          )

          for (col_name in names(reject_cols)) {
            row[[col_name]] <- reject_cols[[col_name]]
          }

          results[[length(results) + 1]] <- row
        }
      }
    }
  }

  if (length(results) == 0) return(NULL)
  dplyr::bind_rows(results)
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
#'   - forecasts: Tibble with OOS forecasts at origin level (NULL if save_forecasts = FALSE)
#' @export
#' @importFrom purrr map_dfr
#' @importFrom dplyr bind_rows mutate filter distinct
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

  # Validate that all series exist in targets_list
  available_series <- names(dataset$targets_list)
  missing_series <- setdiff(series_vec, available_series)

  if (length(missing_series) > 0) {
    log_warn(sprintf(
      "The following series were requested but not found in targets_list (likely removed during data cleaning): %s",
      paste(missing_series, collapse = ", ")
    ), config)
    log_warn(sprintf(
      "Skipping %d series. Continuing with %d available series.",
      length(missing_series), length(series_vec) - length(missing_series)
    ), config)

    # Filter to only valid series
    series_vec <- intersect(series_vec, available_series)
  }

  if (length(series_vec) == 0) {
    stop("No valid series to evaluate after filtering. Check your data cleaning settings.")
  }

  # Get test parameters from config
  test_types <- if (!is.null(config$test_types)) config$test_types else c("DM", "CW")
  test_alpha <- if (!is.null(config$test_alpha)) config$test_alpha else c(0.10, 0.05, 0.01)
  do_tests <- if (!is.null(config$do_tests)) config$do_tests else TRUE

  # Get forecast persistence flag

  do_save_forecasts <- isTRUE(config$save_forecasts)

  # Calculate total iterations for progress tracking
  total_iterations <- length(series_vec) * length(config$horizons) * length(config$factor_methods)
  current_iteration <- 0

  # Initialize result storage
  rmse_list <- list()
  tests_list <- list()
  forecasts_list <- if (do_save_forecasts) list() else NULL

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

        # ========================================
        # EXTRACT FORECASTS FOR MCS PERSISTENCE
        # ========================================
        if (do_save_forecasts) {
          forecast_chunk <- extract_forecasts_from_res(
            res = res,
            dates = dataset$dates,
            series_id = s,
            h = h,
            factor_method = fm,
            eval_start = config$eval_start,
            eval_end = config$eval_end,
            k_max_factor = k_max_factor,
            run_id = config$run_id
          )
          if (!is.null(forecast_chunk)) {
            forecasts_list[[length(forecasts_list) + 1]] <- forecast_chunk
          }
        }
      }
    }
  }

  # Print newline after progress is complete
  cat("\n")

  # Combine all results
  rmse_results <- dplyr::bind_rows(rmse_list)
  tests_results <- if (length(tests_list) > 0) dplyr::bind_rows(tests_list) else NULL

  # Process forecasts: bind and deduplicate AR rows
  forecasts_all <- NULL
  if (do_save_forecasts && length(forecasts_list) > 0) {
    forecasts_all <- dplyr::bind_rows(forecasts_list)

    # Deduplicate AR rows (AR is computed identically for each factor_method)
    # AR rows have factor_method = NA, so they're duplicated across PCA/PLS iterations
    ar_rows <- forecasts_all %>%
      dplyr::filter(model_class == "AR")
    non_ar_rows <- forecasts_all %>%
      dplyr::filter(model_class != "AR")

    # Keep only unique AR rows (by run_id, series_id, h, scheme, origin_index)
    ar_deduped <- ar_rows %>%
      dplyr::distinct(run_id, series_id, h, scheme, origin_index, .keep_all = TRUE)

    forecasts_all <- dplyr::bind_rows(ar_deduped, non_ar_rows)

    log_info(sprintf("Forecasts collected: %d rows (after AR deduplication)",
                     nrow(forecasts_all)), config)
  }

  log_info(sprintf("Unified evaluation complete: RMSE=%d rows, Tests=%d rows",
                   nrow(rmse_results),
                   if (!is.null(tests_results)) nrow(tests_results) else 0), config)

  list(
    rmse_results = rmse_results,
    tests_results = tests_results,
    forecasts = forecasts_all
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

# ============================================================================
# MCS Helper Functions
# ============================================================================

#' Build loss matrix for MCS from saved forecasts
#'
#' Loads persisted forecasts and constructs a loss matrix L suitable for
#' Model Confidence Set (MCS) analysis. The matrix has dimensions (T_oos × M)
#' where T_oos is the number of forecast origins and M is the number of methods.
#'
#' @param forecasts_path Path to forecasts_long.csv, forecasts_long.parquet,
#'   or the forecasts/ directory (partitioned parquet)
#' @param series_id Character: which series to filter
#' @param h Integer: which horizon to filter
#' @param loss_fn Function(y_true, y_hat) -> loss (default: squared error)
#' @param methods_include Character vector of method_ids to include (NULL = all)
#' @return List with:
#'   - L: Matrix (T_oos × M) of losses, balanced panel (complete cases only)
#'   - methods: Character vector of method names (column names of L)
#'   - origin_dates: Date vector of forecast origins corresponding to rows
#'   - n_dropped: Number of rows dropped due to NA in at least one method
#'   - n_total: Total number of rows before dropping
#' @export
#' @importFrom dplyr filter mutate select
#' @importFrom tidyr pivot_wider
#' @importFrom stats complete.cases
build_mcs_loss_matrix <- function(forecasts_path, series_id, h,
                                   loss_fn = function(y, yhat) (y - yhat)^2,
                                   methods_include = NULL) {

  # Load forecasts based on path type
  if (dir.exists(forecasts_path)) {
    # Partitioned parquet directory
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("Package 'arrow' is required to read partitioned parquet. Install with: install.packages('arrow')")
    }
    forecasts <- arrow::open_dataset(forecasts_path) %>%
      dplyr::filter(series_id == !!series_id, h == !!h) %>%
      dplyr::collect()
  } else if (grepl("\\.parquet$", forecasts_path, ignore.case = TRUE)) {
    # Single parquet file
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("Package 'arrow' is required to read parquet files. Install with: install.packages('arrow')")
    }
    forecasts <- arrow::read_parquet(forecasts_path) %>%
      dplyr::filter(series_id == !!series_id, h == !!h)
  } else {
    # CSV file
    forecasts <- readr::read_csv(forecasts_path, show_col_types = FALSE) %>%
      dplyr::filter(series_id == !!series_id, h == !!h)
  }

  if (nrow(forecasts) == 0) {
    stop(sprintf("No forecasts found for series_id='%s', h=%d", series_id, h))
  }

  # Filter to specific methods if requested
  if (!is.null(methods_include)) {
    forecasts <- forecasts %>%
      dplyr::filter(method_id %in% methods_include)
  }

  # Filter to rows with valid y_true and y_hat, then compute loss

  df <- forecasts %>%
    dplyr::filter(!is.na(y_true) & !is.na(y_hat)) %>%
    dplyr::mutate(loss = loss_fn(y_true, y_hat)) %>%
    dplyr::select(origin_index, origin_date, method_id, loss)

  # Pivot to wide format: rows = origin_index, cols = method_id
  L_wide <- df %>%
    tidyr::pivot_wider(
      id_cols = c(origin_index, origin_date),
      names_from = method_id,
      values_from = loss
    ) %>%
    dplyr::arrange(origin_index)

  # Extract method columns (everything except origin_index and origin_date)
  method_cols <- setdiff(names(L_wide), c("origin_index", "origin_date"))

  # Apply complete cases filter (MCS requires balanced panel)
  L_matrix_raw <- as.matrix(L_wide[, method_cols, drop = FALSE])
  complete_mask <- stats::complete.cases(L_matrix_raw)

  n_total <- nrow(L_wide)
  n_dropped <- sum(!complete_mask)

  if (n_dropped / n_total > 0.05) {
    warning(sprintf(
      "Dropped %d/%d rows (%.1f%%) due to NA in at least one method. MCS requires balanced panel.",
      n_dropped, n_total, 100 * n_dropped / n_total
    ))
  }

  L <- L_matrix_raw[complete_mask, , drop = FALSE]
  origin_dates <- L_wide$origin_date[complete_mask]

  list(
    L = L,
    methods = method_cols,
    origin_dates = origin_dates,
    n_dropped = n_dropped,
    n_total = n_total
  )
}
