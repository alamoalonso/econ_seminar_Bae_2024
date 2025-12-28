#' PLS Factor Extraction
#'
#' Functions for extracting factors via Partial Least Squares.
#'
#' @name factors_pls
NULL

#' Extract PLS factors (temporally-aware version)
#'
#' Fits PLS on a subset of data (X_fit, y_fit) but computes scores for all X_all.
#' This is used to enforce temporal constraints in pseudo-out-of-sample forecasting.
#'
#' @param X_fit Matrix (T_fit x N) of predictors for fitting
#' @param y_fit Numeric vector (length T_fit) of target variable for fitting
#' @param X_all Matrix (T_all x N) of predictors for scoring (T_all >= T_fit)
#' @param k_max Maximum number of PLS components to extract
#' @param config Configuration list
#' @return A list with components:
#'   - F: Matrix (T_all x k_max) of PLS scores for all X_all observations
#'   - model: The plsr model object
#' @export
extract_pls_temporal <- function(X_fit, y_fit, X_all, k_max = 12, config = NULL) {
  # Fit PLS on the valid subset
  pls_result_fit <- extract_pls(X_fit, y_fit, k_max = k_max, config = config)

  # Predict scores for all observations
  df_all <- data.frame(y = rep(NA_real_, nrow(X_all)), X_all)
  K <- min(k_max, ncol(X_fit))
  scores_all <- predict(pls_result_fit$model, newdata = df_all[, -1, drop = FALSE], ncomp = 1:K, type = "scores")

  # Convert to matrix
  if (length(dim(scores_all)) == 3) {
    F_all <- scores_all[, , 1, drop = TRUE]
  } else if (length(dim(scores_all)) == 2) {
    F_all <- scores_all
  } else {
    stop("Unexpected dimensions of PLS scores array.")
  }

  F_all <- as.matrix(F_all)
  colnames(F_all) <- paste0("F", seq_len(ncol(F_all)))

  list(
    F = F_all,
    model = pls_result_fit$model
  )
}

#' Extract PLS factors
#'
#' @param X Matrix (T x N) of predictors
#' @param y Numeric vector (length T) of target variable
#' @param k_max Maximum number of PLS components to extract
#' @param config Configuration list
#' @return A list with components:
#'   - F: Matrix (T x k_max) of PLS scores for all observations
#'   - model: The plsr model object
#' @export
#' @importFrom pls plsr
#' @examples
#' X <- matrix(rnorm(100*50), 100, 50)
#' y <- rnorm(100)
#' pls_result <- extract_pls(X, y, k_max = 5, config = config_us_default())
extract_pls <- function(X, y, k_max = 12, config = NULL) {
  # Input validation
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }

  if (length(y) != nrow(X)) {
    stop(sprintf("Dimension mismatch: y has length %d but X has %d rows", length(y), nrow(X)))
  }

  # Build data frame for plsr
  df <- data.frame(y = y, X)

  # Fit only on non-NA observations of y
  df_fit <- df[!is.na(df$y), , drop = FALSE]

  if (nrow(df_fit) <= 2) {
    stop("Too few observations for PLS at this origin.")
  }

  # Key sanity checks for debugging
  if (isTRUE(config$debug)) {
    y_sd <- sd(df_fit$y, na.rm = TRUE)
    log_debug(sprintf(
      "[PLS] nrow(df_fit)=%d, ncol(X)=%d, k_max=%d, sd(y)=%.6f, anyNA(X_fit)=%s",
      nrow(df_fit), ncol(X), k_max, y_sd, anyNA(df_fit[, -1, drop=FALSE])
    ), config)
  }

  # Fit PLS model
  pls_fit <- pls::plsr(
    y ~ .,
    data   = df_fit,
    ncomp  = min(k_max, ncol(X)),
    center = isTRUE(config$pls_center),
    scale  = isTRUE(config$pls_scale),
    method = config$pls_method
  )

  # Diagnostic logging
  if (isTRUE(config$debug)) {
    log_debug(sprintf(
      "[PLS] pls_fit$ncomp=%s | dim(pls_fit$scores)=%s | dim(pls_fit$loading.weights)=%s",
      paste0(pls_fit$ncomp, collapse=""),
      paste(dim(pls_fit$scores), collapse="x"),
      paste(dim(pls_fit$loading.weights), collapse="x")
    ), config)
  }

  # Predict scores for all observations (including those with NA y)
  X_all <- df[, setdiff(names(df), "y"), drop = FALSE]
  K <- min(k_max, ncol(X))
  scores_arr <- predict(pls_fit, newdata = X_all, ncomp = 1:K, type = "scores")

  if (isTRUE(config$debug)) {
    log_debug(sprintf("[PLS] requested ncomp = %s", paste(1:K, collapse=",")), config)
    log_debug(sprintf("[PLS] dim(scores_arr) = %s", paste(dim(scores_arr), collapse="x")), config)
  }

  # Convert to matrix deterministically
  # Handle different array structures returned by predict.mvr
  if (length(dim(scores_arr)) == 3) {
    F_hat <- scores_arr[, , 1, drop = TRUE]
  } else if (length(dim(scores_arr)) == 2) {
    F_hat <- scores_arr
  } else {
    stop("Unexpected dimensions of PLS scores array.")
  }

  F_hat <- as.matrix(F_hat)

  if (isTRUE(config$debug)) {
    log_debug(sprintf(
      "[PLS] dim(F_hat)=%s | ncol(F_hat)=%d | first-row=%s",
      paste(dim(F_hat), collapse="x"),
      ncol(F_hat),
      paste(round(F_hat[1, 1:min(5, ncol(F_hat))], 4), collapse=",")
    ), config)
  }

  # Assign column names
  colnames(F_hat) <- paste0("F", seq_len(ncol(F_hat)))

  # Validate output dimensions
  if (nrow(F_hat) != nrow(X)) {
    stop(sprintf("PLS output dimension mismatch: expected %d rows, got %d", nrow(X), nrow(F_hat)))
  }
  if (ncol(F_hat) != K) {
    stop(sprintf("PLS output dimension mismatch: expected %d columns, got %d", K, ncol(F_hat)))
  }

  list(
    F = F_hat,
    model = pls_fit
  )
}

#' Extract factors at a specific time origin
#'
#' Wrapper that extracts both PCA and PLS factors up to a given time origin.
#'
#' @param panel_final data.frame with date + balanced predictors
#' @param targets_list Nested list of target series
#' @param target_name Character string naming the target series
#' @param h Forecast horizon
#' @param t_origin Date object representing the forecast origin
#' @param k_max_pca Maximum number of PCA components
#' @param k_max_pls Maximum number of PLS components
#' @param config Configuration list
#' @return A list with components pca and pls, each containing F, idx, and other metadata
#' @export
extract_factors_at_origin <- function(panel_final,
                                      targets_list,
                                      target_name,
                                      h,
                                      t_origin,
                                      k_max_pca = 12,
                                      k_max_pls = 12,
                                      config = NULL) {
  # Window up to t_origin
  idx_win <- which(panel_final$date <= t_origin)

  if (length(idx_win) == 0) {
    stop(sprintf("No observations available up to t_origin %s", as.character(t_origin)))
  }

  t_idx <- max(idx_win)  # Current time index

  X_win <- panel_final[idx_win, -1, drop = FALSE]

  # --- PCA ---
  # PCA is unsupervised, so can use all X up to t_origin
  # Dispatch to paper-compliant or library implementation based on config
  if (isTRUE(config$use_paper_pca)) {
    if (isTRUE(config$debug)) {
      log_debug("[extract_factors_at_origin] Using paper-compliant PCA implementation", config)
    }
    pca_out <- extract_pca_paper(X_win, k_max = k_max_pca, config = config, use_paper_method = TRUE)
  } else {
    if (isTRUE(config$debug)) {
      log_debug("[extract_factors_at_origin] Using library PCA (prcomp)", config)
    }
    pca_out <- extract_pca(X_win, k_max = k_max_pca, config = config)
  }

  # --- PLS ---
  # PLS is supervised: must respect temporal constraints
  # At time t_idx, we can only observe y_h[t] where t + h <= t_idx
  # This means we can only use training data for t <= t_idx - h
  y_full <- targets_list[[target_name]][[paste0("h", h)]]

  # CRITICAL FIX: Filter to only use y values that are observed at t_origin
  # For y_h[t] = x[t+h] to be known at time t_idx, we need t+h <= t_idx
  valid_pls_idx <- idx_win[idx_win + h <= t_idx]

  if (isTRUE(config$debug)) {
    log_debug(sprintf(
      "[extract_factors_at_origin] PLS temporal filtering: t_idx=%d, h=%d, valid range=1:%d (was 1:%d)",
      t_idx, h, max(valid_pls_idx), t_idx
    ), config)
  }

  if (length(valid_pls_idx) < 2) {
    # Not enough observations for PLS estimation
    log_warn(sprintf(
      "[extract_factors_at_origin] Insufficient observations for PLS at t_idx=%d, h=%d (need at least 2)",
      t_idx, h
    ), config)
    # Return factors with NAs
    pls_out <- list(
      F = matrix(NA_real_, nrow = length(idx_win), ncol = k_max_pls),
      model = NULL
    )
  } else {
    # Extract valid training data for PLS
    # PLS will be fit using (X_pls_fit, y_pls_fit) but scores computed for all X_win
    X_pls_fit <- as.matrix(panel_final[valid_pls_idx, -1, drop = FALSE])
    y_pls_fit <- y_full[valid_pls_idx]

    # Dispatch to paper-compliant or library implementation based on config
    if (isTRUE(config$use_paper_pls)) {
      if (isTRUE(config$debug)) {
        log_debug("[extract_factors_at_origin] Using paper-compliant PLS implementation", config)
      }
      # Paper PLS: fit on valid data, compute scores for all X_win
      # Need to pass both fitting data and full data
      pls_out <- extract_pls_paper_temporal(
        X_fit = X_pls_fit,
        y_fit = y_pls_fit,
        X_all = X_win,
        k_max = k_max_pls,
        config = config,
        use_paper_method = TRUE
      )
    } else {
      if (isTRUE(config$debug)) {
        log_debug("[extract_factors_at_origin] Using library PLS (pls::plsr)", config)
      }
      # Library PLS: fit on valid data, predict scores for all X_win
      pls_out <- extract_pls_temporal(
        X_fit = X_pls_fit,
        y_fit = y_pls_fit,
        X_all = X_win,
        k_max = k_max_pls,
        config = config
      )
    }
  }

  list(
    pca = list(
      F        = pca_out$F,        # length(idx_win) x k_pca
      loadings = pca_out$loadings,
      idx      = idx_win
    ),
    pls = list(
      F     = pls_out$F,           # length(idx_win) x k_pls
      model = pls_out$model,
      idx   = idx_win
    )
  )
}
