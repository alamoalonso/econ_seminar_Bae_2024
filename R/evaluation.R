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
