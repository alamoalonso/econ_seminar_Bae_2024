#' Forecasting Model Runner
#'
#' High-level functions to run forecasts for a series across all time origins.
#'
#' @name forecasting_models
NULL

#' Run forecasts for a single series and horizon
#'
#' Performs rolling/recursive forecasting for one target series at one horizon.
#' Uses a single factor method (PCA or PLS) and generates predictions for all models.
#'
#' @param series_name Character string: name of the target series
#' @param h Integer: forecast horizon
#' @param dates Vector of Date objects
#' @param panel_final data.frame with balanced predictors
#' @param panel_std1 data.frame with all standardized predictors (for y_base)
#' @param targets_list Nested list of target series
#' @param factor_method Character: "PCA" or "PLS"
#' @param k_max_pca Integer: maximum number of PCA components
#' @param k_max_pls Integer: maximum number of PLS components
#' @param config Configuration list
#' @return A list with components:
#'   - truth: vector of true values
#'   - DI_rec, DI_roll: matrices (T x k_max) of DI predictions
#'   - DIAR_rec, DIAR_roll: matrices (T x k_max) of DIAR predictions
#'   - DLAG_rec, DLAG_roll: matrices (T x min(4, k_max)) of DIAR-LAG predictions
#'   - AR_rec, AR_roll: vectors of AR predictions
#' @export
run_forecasts_for_series <- function(series_name,
                                     h,
                                     dates,
                                     panel_final,
                                     panel_std1,
                                     targets_list,
                                     factor_method = c("PCA", "PLS"),
                                     k_max_pca = 12,
                                     k_max_pls = 12,
                                     config = NULL) {

  factor_method <- match.arg(factor_method)

  log_debug(sprintf(
    "Running forecasts: series=%s, h=%d, factor_method=%s",
    series_name, h, factor_method
  ), config)

  y_h    <- targets_list[[series_name]][[paste0("h", h)]]
  y_base <- panel_std1[[series_name]]

  Tn <- length(y_h)

  # Determine k_max based on factor method
  k_max_common <- if (factor_method == "PCA") k_max_pca else k_max_pls

  # Initialize storage
  out <- list(
    truth    = rep(NA_real_, Tn),
    DI_rec   = matrix(NA, Tn, k_max_common),
    DI_roll  = matrix(NA, Tn, k_max_common),
    DIAR_rec = matrix(NA, Tn, k_max_common),
    DIAR_roll= matrix(NA, Tn, k_max_common),
    DLAG_rec = matrix(NA, Tn, min(4, k_max_common)),
    DLAG_roll= matrix(NA, Tn, min(4, k_max_common)),
    AR_rec   = rep(NA_real_, Tn),
    AR_roll  = rep(NA_real_, Tn)
  )

  # Forecasting loop
  first_idx <- config$first_forecast_idx
  if (is.null(first_idx)) first_idx <- 60

  for (t_idx in first_idx:(Tn - h)) {

    # Define true value
    out$truth[t_idx] <- y_h[t_idx]
    if (is.na(out$truth[t_idx])) next

    t_origin <- dates[t_idx]

    # Extract factors at this origin
    fcts <- extract_factors_at_origin(
      panel_final  = panel_final,
      targets_list = targets_list,
      target_name  = series_name,
      h            = h,
      t_origin     = t_origin,
      k_max_pca    = k_max_pca,
      k_max_pls    = k_max_pls,
      config       = config
    )

    # Select factor matrix based on method
    if (factor_method == "PCA") {
      F_used <- fcts$pca$F
    } else {
      F_used <- fcts$pls$F
    }

    # Debug logging at trace origins
    if (should_trace(t_idx, config)) {
      log_debug(sprintf(
        "[RUN] t_idx=%d date=%s | dim(F_used)=%s | k_max_common=%d | k_max_used=%d",
        t_idx, as.character(t_origin),
        paste(dim(F_used), collapse="x"),
        k_max_common, min(k_max_common, ncol(F_used))
      ), config)
    }

    # Safety: ensure k_max does not exceed available components
    k_max_used <- min(k_max_common, ncol(F_used))

    # Define sample indices
    idx_rec  <- get_recursive_idx(t_idx)
    idx_roll <- get_rolling_idx(dates, t_idx)

    # -------- AR (BIC) --------
    ar_fit_rec  <- fit_AR_BIC(y_h, y_base, idx_rec, max_p = config$max_p_ar)
    ar_fit_roll <- fit_AR_BIC(y_h, y_base, idx_roll, max_p = config$max_p_ar)

    out$AR_rec[t_idx]  <- predict_AR_at_origin(ar_fit_rec,  y_base, t_idx)
    out$AR_roll[t_idx] <- predict_AR_at_origin(ar_fit_roll, y_base, t_idx)

    # -------- DI --------
    # TEMPORAL FIX: Pass h and t_idx to enforce t + h <= t_idx constraint
    for (k in 1:k_max_used) {
      di_fit_rec  <- fit_DI(y_h, F_used, k, idx_rec, h = h, t_current = t_idx)
      di_fit_roll <- fit_DI(y_h, F_used, k, idx_roll, h = h, t_current = t_idx)

      if (!is.null(di_fit_rec))
        out$DI_rec[t_idx, k] <- predict_DI_at_origin(di_fit_rec, F_used, t_idx, k)
      if (!is.null(di_fit_roll))
        out$DI_roll[t_idx, k] <- predict_DI_at_origin(di_fit_roll, F_used, t_idx, k)
    }

    # -------- DIAR --------
    # TEMPORAL FIX: Pass h and t_idx to enforce t + h <= t_idx constraint
    for (k in 1:k_max_used) {
      diar_fit_rec  <- fit_DIAR_BIC(y_h, y_base, F_used, k, idx_rec, max_p = config$max_p_ar, h = h, t_current = t_idx)
      diar_fit_roll <- fit_DIAR_BIC(y_h, y_base, F_used, k, idx_roll, max_p = config$max_p_ar, h = h, t_current = t_idx)

      if (!is.null(diar_fit_rec))
        out$DIAR_rec[t_idx, k] <- predict_DIAR_at_origin(diar_fit_rec, F_used, y_base, t_idx, k)
      if (!is.null(diar_fit_roll))
        out$DIAR_roll[t_idx, k] <- predict_DIAR_at_origin(diar_fit_roll, F_used, y_base, t_idx, k)
    }

    # -------- DIAR-LAG (k <= 4) --------
    # TEMPORAL FIX: Pass h and t_idx to enforce t + h <= t_idx constraint
    for (k in 1:min(4, k_max_used)) {
      dlag_fit_rec  <- fit_DIARLAG_BIC(y_h, y_base, F_used, k, idx_rec,
                                       max_p = config$max_p_ar, max_m = config$max_m_diarlag, h = h, t_current = t_idx)
      dlag_fit_roll <- fit_DIARLAG_BIC(y_h, y_base, F_used, k, idx_roll,
                                       max_p = config$max_p_ar, max_m = config$max_m_diarlag, h = h, t_current = t_idx)

      if (!is.null(dlag_fit_rec))
        out$DLAG_rec[t_idx, k] <- predict_DIARLAG_at_origin(dlag_fit_rec, F_used, y_base, t_idx, k)
      if (!is.null(dlag_fit_roll))
        out$DLAG_roll[t_idx, k] <- predict_DIARLAG_at_origin(dlag_fit_roll, F_used, y_base, t_idx, k)
    }
  }

  log_debug(sprintf(
    "Completed forecasts: series=%s, h=%d, factor_method=%s",
    series_name, h, factor_method
  ), config)

  out
}

#' Run forecast series with explicit configuration (wrapper)
#'
#' This is the main entry point for running forecasts for a single series.
#'
#' @param series_id Character string: name of the series to forecast
#' @param dataset List returned from preprocess_dataset()
#' @param config Configuration list
#' @return List of forecast results for all horizons and factor methods
#' @export
run_forecast_series <- function(series_id, dataset, config) {
  log_debug(sprintf("Running forecast series: %s", series_id), config)

  results <- list()

  for (h in config$horizons) {
    for (fm in config$factor_methods) {
      key <- sprintf("h%d_%s", h, fm)

      results[[key]] <- run_forecasts_for_series(
        series_name   = series_id,
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
    }
  }

  results
}
