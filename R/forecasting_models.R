#' Forecasting Model Runner
#'
#' High-level functions to run forecasts for a series across all time origins.
#' Supports both grid (fixed k=1..k_max) and dynamic (k_hat selected per origin)
#' factor specifications in the same workflow run.
#'
#' @name forecasting_models
NULL

# Helper: null-coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

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
    idx_roll <- get_rolling_idx(dates, t_idx, config)

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

# ============================================================================
# Factor Spec-Based Forecasting (Mixed Grid + Dynamic)
# ============================================================================

#' Run forecasts for a single series, horizon, and factor_spec
#'
#' Performs rolling/recursive forecasting using a specific factor specification.
#' For grid specs: evaluates k=1..k_max
#' For dynamic specs: selects k_hat(t) at each origin using the specified rule
#'
#' @param series_name Character string: name of the target series
#' @param h Integer: forecast horizon
#' @param dates Vector of Date objects
#' @param panel_final data.frame with balanced predictors
#' @param panel_std1 data.frame with all standardized predictors (for y_base)
#' @param targets_list Nested list of target series
#' @param factor_spec A factor_spec object (from create_factor_spec or get_factor_specs)
#' @param config Configuration list
#' @return A list with components organized by the spec:
#'   - spec_id: The factor_spec id
#'   - k_mode: "grid" or "dynamic"
#'   - k_rule: NULL for grid, rule name for dynamic
#'   - truth: vector of true values
#'   - k_hat: vector of selected k per origin (for dynamic; NA for grid)
#'   - DI_rec, DI_roll: For grid: matrices (T x k_max); For dynamic: vectors
#'   - DIAR_rec, DIAR_roll: Same structure
#'   - DLAG_rec, DLAG_roll: Same structure (k <= 4)
#'   - AR_rec, AR_roll: vectors of AR predictions
#' @export
run_forecasts_for_spec <- function(series_name,
                                    h,
                                    dates,
                                    panel_final,
                                    panel_std1,
                                    targets_list,
                                    factor_spec,
                                    config = NULL) {

  spec_id <- factor_spec$id
  factor_method <- factor_spec$factor_method
  k_mode <- factor_spec$k_mode
  k_rule <- factor_spec$k_rule
  k_max <- factor_spec$k_max

  log_debug(sprintf(
    "Running forecasts: series=%s, h=%d, spec=%s (%s, %s)",
    series_name, h, spec_id, factor_method, k_mode
  ), config)

  y_h <- targets_list[[series_name]][[paste0("h", h)]]
  y_base <- panel_std1[[series_name]]

  Tn <- length(y_h)

  # Initialize storage based on k_mode
  if (k_mode == "grid") {
    out <- list(
      spec_id = spec_id,
      factor_method = factor_method,
      k_mode = k_mode,
      k_rule = k_rule,
      k_max = k_max,
      truth = rep(NA_real_, Tn),
      k_hat = rep(NA_integer_, Tn),  # Not used for grid, but included for consistency
      DI_rec = matrix(NA, Tn, k_max),
      DI_roll = matrix(NA, Tn, k_max),
      DIAR_rec = matrix(NA, Tn, k_max),
      DIAR_roll = matrix(NA, Tn, k_max),
      DLAG_rec = matrix(NA, Tn, min(4, k_max)),
      DLAG_roll = matrix(NA, Tn, min(4, k_max)),
      AR_rec = rep(NA_real_, Tn),
      AR_roll = rep(NA_real_, Tn)
    )
  } else {
    # Dynamic: single forecast per origin (k = k_hat)
    out <- list(
      spec_id = spec_id,
      factor_method = factor_method,
      k_mode = k_mode,
      k_rule = k_rule,
      k_max = k_max,
      truth = rep(NA_real_, Tn),
      k_hat = rep(NA_integer_, Tn),       # Selected k at each origin (for DI/DIAR)
      k_hat_dlag = rep(NA_integer_, Tn),  # Selected k for DIAR-LAG (k_max = 4)
      DI_rec = rep(NA_real_, Tn),
      DI_roll = rep(NA_real_, Tn),
      DIAR_rec = rep(NA_real_, Tn),
      DIAR_roll = rep(NA_real_, Tn),
      DLAG_rec = rep(NA_real_, Tn),
      DLAG_roll = rep(NA_real_, Tn),
      AR_rec = rep(NA_real_, Tn),
      AR_roll = rep(NA_real_, Tn)
    )
  }

  # Forecasting loop
  first_idx <- config$first_forecast_idx %||% 60

  for (t_idx in first_idx:(Tn - h)) {

    # Define true value
    out$truth[t_idx] <- y_h[t_idx]
    if (is.na(out$truth[t_idx])) next

    t_origin <- dates[t_idx]

    # --- Factor Extraction ---
    # Use the standard extraction which respects temporal constraints
    fcts <- extract_factors_at_origin(
      panel_final = panel_final,
      targets_list = targets_list,
      target_name = series_name,
      h = h,
      t_origin = t_origin,
      k_max_pca = if (factor_method == "PCA") k_max else config$k_max_pca,
      k_max_pls = if (factor_method %in% c("PLS", "1-PLS")) k_max else config$k_max_pls,
      config = config
    )

    # Select factor matrix based on method
    if (factor_method == "PCA") {
      F_used <- fcts$pca$F
    } else {
      F_used <- fcts$pls$F
    }

    # Safety: ensure k_max does not exceed available components
    k_max_used <- min(k_max, ncol(F_used))

    # Define sample indices
    idx_rec <- get_recursive_idx(t_idx)
    idx_roll <- get_rolling_idx(dates, t_idx, config)

    # --- AR (BIC) --- same for all specs
    ar_fit_rec <- fit_AR_BIC(y_h, y_base, idx_rec, max_p = config$max_p_ar)
    ar_fit_roll <- fit_AR_BIC(y_h, y_base, idx_roll, max_p = config$max_p_ar)

    out$AR_rec[t_idx] <- predict_AR_at_origin(ar_fit_rec, y_base, t_idx)
    out$AR_roll[t_idx] <- predict_AR_at_origin(ar_fit_roll, y_base, t_idx)

    # --- Determine k value(s) to evaluate ---
    if (k_mode == "grid") {
      k_values <- 1:k_max_used
    } else {
      # Dynamic: compute k_hat for DI/DIAR (using full k_max)
      k_hat <- compute_k_hat_at_origin(
        panel_final = panel_final,
        t_idx = t_idx,
        factor_method = factor_method,
        k_rule = k_rule,
        k_max = k_max_used,
        h = h,
        config = config,
        context = list(series_id = series_name, origin_index = t_idx, origin_date = t_origin)
      )
      out$k_hat[t_idx] <- k_hat
      k_values <- k_hat

      # Dynamic: compute k_hat_dlag for DIAR-LAG (using k_max = 4, per Bae 2024)
      k_max_dlag <- min(4, k_max_used)
      k_hat_dlag <- compute_k_hat_at_origin(
        panel_final = panel_final,
        t_idx = t_idx,
        factor_method = factor_method,
        k_rule = k_rule,
        k_max = k_max_dlag,
        h = h,
        config = config,
        context = list(series_id = series_name, origin_index = t_idx, origin_date = t_origin, model = "DIAR-LAG")
      )
      out$k_hat_dlag[t_idx] <- k_hat_dlag
    }

    # --- Generate forecasts for each k in k_values ---
    for (k in k_values) {
      # DI
      di_fit_rec <- fit_DI(y_h, F_used, k, idx_rec, h = h, t_current = t_idx)
      di_fit_roll <- fit_DI(y_h, F_used, k, idx_roll, h = h, t_current = t_idx)

      if (k_mode == "grid") {
        if (!is.null(di_fit_rec))
          out$DI_rec[t_idx, k] <- predict_DI_at_origin(di_fit_rec, F_used, t_idx, k)
        if (!is.null(di_fit_roll))
          out$DI_roll[t_idx, k] <- predict_DI_at_origin(di_fit_roll, F_used, t_idx, k)
      } else {
        if (!is.null(di_fit_rec))
          out$DI_rec[t_idx] <- predict_DI_at_origin(di_fit_rec, F_used, t_idx, k)
        if (!is.null(di_fit_roll))
          out$DI_roll[t_idx] <- predict_DI_at_origin(di_fit_roll, F_used, t_idx, k)
      }

      # DIAR
      diar_fit_rec <- fit_DIAR_BIC(y_h, y_base, F_used, k, idx_rec, max_p = config$max_p_ar, h = h, t_current = t_idx)
      diar_fit_roll <- fit_DIAR_BIC(y_h, y_base, F_used, k, idx_roll, max_p = config$max_p_ar, h = h, t_current = t_idx)

      if (k_mode == "grid") {
        if (!is.null(diar_fit_rec))
          out$DIAR_rec[t_idx, k] <- predict_DIAR_at_origin(diar_fit_rec, F_used, y_base, t_idx, k)
        if (!is.null(diar_fit_roll))
          out$DIAR_roll[t_idx, k] <- predict_DIAR_at_origin(diar_fit_roll, F_used, y_base, t_idx, k)
      } else {
        if (!is.null(diar_fit_rec))
          out$DIAR_rec[t_idx] <- predict_DIAR_at_origin(diar_fit_rec, F_used, y_base, t_idx, k)
        if (!is.null(diar_fit_roll))
          out$DIAR_roll[t_idx] <- predict_DIAR_at_origin(diar_fit_roll, F_used, y_base, t_idx, k)
      }

      # DIAR-LAG (only for k <= 4) - grid mode only within this loop
      if (k_mode == "grid" && k <= min(4, k_max_used)) {
        dlag_fit_rec <- fit_DIARLAG_BIC(y_h, y_base, F_used, k, idx_rec,
                                         max_p = config$max_p_ar, max_m = config$max_m_diarlag, h = h, t_current = t_idx)
        dlag_fit_roll <- fit_DIARLAG_BIC(y_h, y_base, F_used, k, idx_roll,
                                          max_p = config$max_p_ar, max_m = config$max_m_diarlag, h = h, t_current = t_idx)

        if (!is.null(dlag_fit_rec))
          out$DLAG_rec[t_idx, k] <- predict_DIARLAG_at_origin(dlag_fit_rec, F_used, y_base, t_idx, k)
        if (!is.null(dlag_fit_roll))
          out$DLAG_roll[t_idx, k] <- predict_DIARLAG_at_origin(dlag_fit_roll, F_used, y_base, t_idx, k)
      }
    }

    # DIAR-LAG for dynamic mode: use k_hat_dlag (computed with k_max = 4)
    if (k_mode == "dynamic") {
      k_dlag <- out$k_hat_dlag[t_idx]
      if (!is.na(k_dlag) && k_dlag >= 1) {
        dlag_fit_rec <- fit_DIARLAG_BIC(y_h, y_base, F_used, k_dlag, idx_rec,
                                         max_p = config$max_p_ar, max_m = config$max_m_diarlag, h = h, t_current = t_idx)
        dlag_fit_roll <- fit_DIARLAG_BIC(y_h, y_base, F_used, k_dlag, idx_roll,
                                          max_p = config$max_p_ar, max_m = config$max_m_diarlag, h = h, t_current = t_idx)

        if (!is.null(dlag_fit_rec))
          out$DLAG_rec[t_idx] <- predict_DIARLAG_at_origin(dlag_fit_rec, F_used, y_base, t_idx, k_dlag)
        if (!is.null(dlag_fit_roll))
          out$DLAG_roll[t_idx] <- predict_DIARLAG_at_origin(dlag_fit_roll, F_used, y_base, t_idx, k_dlag)
      }
    }
  }

  log_debug(sprintf(
    "Completed forecasts: series=%s, h=%d, spec=%s",
    series_name, h, spec_id
  ), config)

  out
}


#' Compute k_hat at a single origin using dynamic rule
#'
#' @param panel_final data.frame with balanced predictors
#' @param t_idx Current time index (origin)
#' @param factor_method Factor extraction method
#' @param k_rule Decision rule: "bn_bic" or "onatski"
#' @param k_max Maximum k allowed
#' @param h Forecast horizon
#' @param config Configuration list
#' @param context Contextual info for logging
#' @return Integer k_hat
#' @keywords internal
compute_k_hat_at_origin <- function(panel_final, t_idx, factor_method, k_rule, k_max, h, config, context = NULL) {

  # Extract training X up to origin
  idx_win <- 1:t_idx
  X_train <- as.matrix(panel_final[idx_win, -1, drop = FALSE])

  T_obs <- nrow(X_train)
  N <- ncol(X_train)

  # Get settings
  settings <- config$k_selection_settings %||% list()
  min_k <- settings$min_k %||% 1L

  if (k_rule == "bn_bic") {
    # BN-BIC for PCA
    # Apply PCA preprocessing
    pca_center <- isTRUE(config$pca_center)
    pca_scale <- isTRUE(config$pca_scale)
    X_preprocessed <- scale(X_train, center = pca_center, scale = pca_scale)

    sigma_sq_rule <- settings$bn_bic_sigma_sq %||% "v_kmax"
    result <- compute_bn_bic_k(X_preprocessed, k_max, sigma_sq_rule, min_k, config)

    if (!is.null(result$warning) && !is.null(context)) {
      log_warn(sprintf(
        "[BN-BIC] %s | series=%s, origin=%d, date=%s",
        result$warning, context$series_id, t_idx, as.character(context$origin_date)
      ), config)
    }

    return(result$k_hat)

  } else if (k_rule == "onatski") {
    # Onatski for PLS
    # Apply PLS preprocessing (same as used for PLS extraction)
    pls_center <- isTRUE(config$pls_center)
    pls_scale <- isTRUE(config$pls_scale)
    X_preprocessed <- scale(X_train, center = pls_center, scale = pls_scale)

    r_max <- settings$onatski_r_max %||% 12L
    result <- compute_onatski_k(X_preprocessed, r_max, min_k, k_max, config)

    if (!is.null(result$warning) && !is.null(context)) {
      log_warn(sprintf(
        "[Onatski] %s | series=%s, origin=%d, date=%s",
        result$warning, context$series_id, t_idx, as.character(context$origin_date)
      ), config)
    }

    return(result$k_hat)

  } else {
    stop(sprintf("Unknown k_rule: '%s'", k_rule))
  }
}


#' Run forecasts for all factor_specs for a series and horizon
#'
#' Iterates over all factor_specs and runs forecasts for each.
#' This is the main entry point for spec-based forecasting.
#'
#' @param series_name Character string: name of the target series
#' @param h Integer: forecast horizon
#' @param dataset List with panel_final, panel_std1, targets_list, dates
#' @param factor_specs List of factor_spec objects
#' @param config Configuration list
#' @return Named list of forecast results, keyed by spec_id
#' @export
run_forecasts_for_all_specs <- function(series_name, h, dataset, factor_specs, config = NULL) {

  results <- list()

  for (spec in factor_specs) {
    results[[spec$id]] <- run_forecasts_for_spec(
      series_name = series_name,
      h = h,
      dates = dataset$dates,
      panel_final = dataset$panel_final,
      panel_std1 = dataset$panel_std1,
      targets_list = dataset$targets_list,
      factor_spec = spec,
      config = config
    )
  }

  results
}
