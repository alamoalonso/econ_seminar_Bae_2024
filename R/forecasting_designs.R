#' Forecasting Design Functions
#'
#' Functions for building, fitting, and predicting with various forecasting models:
#' - AR: Autoregressive model with BIC-based lag selection
#' - DI: Direct forecast with factors
#' - DIAR: Direct forecast with factors and AR lags
#' - DIAR-LAG: Direct forecast with factors and their lags (k <= 4)
#'
#' @name forecasting_designs
NULL

# ============================================================================
# Utility functions
# ============================================================================

#' Compute BIC
#'
#' @param residuals Numeric vector of model residuals
#' @param n Sample size
#' @param k_param Number of parameters
#' @return BIC value
compute_bic <- function(residuals, n, k_param) {
  rss <- sum(residuals^2)
  sigma2 <- rss / n
  n * log(sigma2) + k_param * log(n)
}

#' Get recursive window indices (1 to t_idx)
#'
#' @param t_idx Current time index
#' @return Integer vector of indices
#' @export
get_recursive_idx <- function(t_idx) {
  1:t_idx
}

#' Get rolling window indices (observation-based window ending at t_idx)
#'
#' Computes a rolling window based on the config$rolling_window_years parameter
#' and the data frequency. Converts the window from calendar years to observations.
#'
#' @param dates Vector of Date objects
#' @param t_idx Current time index
#' @param config Configuration list (must contain frequency and rolling_window_years)
#' @return Integer vector of indices
#' @export
get_rolling_idx <- function(dates, t_idx, config) {
  # Convert rolling_window_years to observations based on frequency
  obs_per_year <- switch(config$frequency,
    "monthly" = 12,
    "quarterly" = 4,
    "yearly" = 1,
    stop(sprintf("Unknown frequency: %s. Must be 'monthly', 'quarterly', or 'yearly'.",
                 config$frequency))
  )

  window_size_obs <- config$rolling_window_years * obs_per_year

  # Use observation-based window (more robust than date arithmetic)
  start_idx <- max(1, t_idx - window_size_obs + 1)
  start_idx:t_idx
}

# ============================================================================
# AR (Autoregressive) model
# ============================================================================

#' Build AR design matrix
#'
#' @param y_h Numeric vector: h-step target series (length T)
#' @param y_base Numeric vector: base series for lags (length T)
#' @param p Integer: number of AR lags
#' @param sample_idx Integer vector: indices to use for estimation
#' @return List with y and X, or NULL if insufficient data
build_design_AR <- function(y_h, y_base, p, sample_idx) {
  max_lag <- p

  # Need t > p for all lags to be defined
  idx <- sample_idx[sample_idx > max_lag & !is.na(y_h[sample_idx])]
  if (length(idx) == 0) return(NULL)

  y <- y_h[idx]

  if (p > 0) {
    X_ylags <- sapply(1:p, function(j) y_base[idx - j])
    colnames(X_ylags) <- paste0("ylag", 1:p)
    X_raw <- X_ylags
  } else {
    # p = 0 => constant only
    X_raw <- NULL
  }

  if (!is.null(X_raw)) {
    ok <- stats::complete.cases(cbind(y, X_raw))
    y <- y[ok]
    X_raw <- X_raw[ok, , drop = FALSE]
  } else {
    ok <- stats::complete.cases(y)
    y <- y[ok]
  }

  if (length(y) == 0) return(NULL)

  if (is.null(X_raw)) {
    X <- matrix(1, nrow = length(y), ncol = 1)
    colnames(X) <- "Intercept"
  } else {
    X <- cbind(Intercept = 1, X_raw)
  }

  list(y = y, X = X)
}

#' Fit AR model with BIC-based lag selection
#'
#' @param y_h Numeric vector: h-step target series
#' @param y_base Numeric vector: base series for lags
#' @param sample_idx Integer vector: indices to use for estimation
#' @param max_p Maximum lag order to consider (default 6)
#' @return List with model, bic, p, n, or NULL if no valid fit
#' @export
fit_AR_BIC <- function(y_h, y_base, sample_idx, max_p = 6) {
  best_bic <- Inf
  best_res <- NULL

  for (p in 0:max_p) {
    des <- build_design_AR(y_h, y_base, p, sample_idx)
    if (is.null(des)) next

    mod <- lm(des$y ~ des$X - 1)
    n <- length(des$y)
    k_param <- ncol(des$X)
    bic <- compute_bic(residuals(mod), n, k_param)

    if (bic < best_bic) {
      best_bic <- bic
      best_res <- list(
        model = mod,
        bic   = bic,
        p     = p,
        n     = n
      )
    }
  }

  best_res
}

#' Build AR predictor row for a given time index
#'
#' @param y_base Numeric vector: base series
#' @param t_index Integer: time index for prediction
#' @param p Integer: number of AR lags
#' @return Matrix (1 x (p+1)) or NULL if insufficient data
build_AR_row_for_origin <- function(y_base, t_index, p) {
  if (p == 0) {
    x_row <- c(Intercept = 1)
  } else {
    if (t_index <= p) return(NULL)
    lags <- sapply(1:p, function(j) y_base[t_index - j])
    if (any(!is.finite(lags))) return(NULL)
    x_row <- c(Intercept = 1, setNames(lags, paste0("ylag", 1:p)))
  }
  matrix(x_row, nrow = 1)
}

#' Predict AR at a given time origin
#'
#' @param fit_obj List returned from fit_AR_BIC
#' @param y_base Numeric vector: base series
#' @param t_index Integer: time index for prediction
#' @return Numeric: predicted value, or NA if prediction fails
#' @export
predict_AR_at_origin <- function(fit_obj, y_base, t_index) {
  if (is.null(fit_obj)) return(NA_real_)
  p <- fit_obj$p
  x_row <- build_AR_row_for_origin(y_base, t_index, p)
  if (is.null(x_row)) return(NA_real_)

  beta <- coef(fit_obj$model)
  as.numeric(x_row %*% beta)
}

# ============================================================================
# DI (Direct forecast with factors)
# ============================================================================

#' Build DI design matrix
#'
#' @param y_h Numeric vector: h-step target series (length T)
#' @param F Matrix: factor matrix (T x Kmax)
#' @param k Integer: number of factors to use
#' @param sample_idx Integer vector: indices to use for estimation
#' @param h Integer: forecast horizon (for temporal constraint checking, optional)
#' @param t_current Integer: current time index (for temporal constraint checking, optional)
#' @return List with y and X, or NULL if insufficient data
build_design_DI <- function(y_h, F, k, sample_idx, h = NULL, t_current = NULL) {
  # TEMPORAL CONSTRAINT: At time t_current, we can only use y_h[t] where t+h <= t_current
  # This means we should only include indices t where t <= t_current - h
  if (!is.null(h) && !is.null(t_current)) {
    # Filter to respect temporal constraint: t + h <= t_current
    sample_idx <- sample_idx[sample_idx + h <= t_current]
  }

  idx <- sample_idx[!is.na(y_h[sample_idx])]
  if (length(idx) == 0) return(NULL)

  y <- y_h[idx]
  X_f <- as.matrix(F[idx, 1:k, drop = FALSE])

  ok <- stats::complete.cases(cbind(y, X_f))
  y <- y[ok]
  X_f <- X_f[ok, , drop = FALSE]

  if (length(y) == 0) return(NULL)

  X <- cbind(Intercept = 1, X_f)
  list(y = y, X = X)
}

#' Fit DI model
#'
#' @param y_h Numeric vector: h-step target series
#' @param F Matrix: factor matrix
#' @param k Integer: number of factors to use
#' @param sample_idx Integer vector: indices to use for estimation
#' @param h Integer: forecast horizon (for temporal constraint)
#' @param t_current Integer: current time index (for temporal constraint)
#' @return List with model, or NULL if no valid fit
#' @export
fit_DI <- function(y_h, F, k, sample_idx, h = NULL, t_current = NULL) {
  des <- build_design_DI(y_h, F, k, sample_idx, h = h, t_current = t_current)
  if (is.null(des)) return(NULL)

  mod <- lm(des$y ~ des$X - 1)
  list(model = mod)
}

#' Predict DI at a given time origin
#'
#' @param fit_obj List returned from fit_DI
#' @param F Matrix: factor matrix
#' @param t_idx Integer: time index for prediction
#' @param k Integer: number of factors to use
#' @return Numeric: predicted value, or NA if prediction fails
#' @export
predict_DI_at_origin <- function(fit_obj, F, t_idx, k) {
  if (is.null(fit_obj)) return(NA_real_)
  beta <- coef(fit_obj$model)

  # Intercept + first k factors
  x <- c(Intercept = 1, as.numeric(F[t_idx, 1:k]))
  if (length(beta) != length(x)) return(NA_real_)

  as.numeric(sum(x * beta))
}

# ============================================================================
# DIAR (Direct forecast with factors and AR lags)
# ============================================================================

#' Build DIAR design matrix
#'
#' @param y_h Numeric vector: h-step target series
#' @param y_base Numeric vector: base series for lags
#' @param F Matrix: factor matrix
#' @param k Integer: number of factors to use
#' @param p Integer: number of AR lags
#' @param sample_idx Integer vector: indices to use for estimation
#' @param h Integer: forecast horizon (for temporal constraint checking, optional)
#' @param t_current Integer: current time index (for temporal constraint checking, optional)
#' @return List with y and X, or NULL if insufficient data
build_design_DIAR <- function(y_h, y_base, F, k, p, sample_idx, h = NULL, t_current = NULL) {
  # TEMPORAL CONSTRAINT: Filter indices to respect t + h <= t_current
  if (!is.null(h) && !is.null(t_current)) {
    sample_idx <- sample_idx[sample_idx + h <= t_current]
  }

  max_lag <- p
  idx <- sample_idx[sample_idx > max_lag & !is.na(y_h[sample_idx])]
  if (length(idx) == 0) return(NULL)

  y <- y_h[idx]
  X_f <- as.matrix(F[idx, 1:k, drop = FALSE])

  if (p > 0) {
    X_ylags <- sapply(1:p, function(j) y_base[idx - j])
    colnames(X_ylags) <- paste0("ylag", 1:p)
    X_raw <- cbind(X_f, X_ylags)
  } else {
    X_raw <- X_f
  }

  ok <- stats::complete.cases(cbind(y, X_raw))
  y <- y[ok]
  X_raw <- X_raw[ok, , drop = FALSE]

  if (length(y) == 0) return(NULL)

  X <- cbind(Intercept = 1, X_raw)
  list(y = y, X = X)
}

#' Fit DIAR model with BIC-based lag selection
#'
#' @param y_h Numeric vector: h-step target series
#' @param y_base Numeric vector: base series for lags
#' @param F Matrix: factor matrix
#' @param k Integer: number of factors to use
#' @param sample_idx Integer vector: indices to use for estimation
#' @param max_p Maximum lag order to consider (default 6)
#' @param h Integer: forecast horizon (for temporal constraint)
#' @param t_current Integer: current time index (for temporal constraint)
#' @return List with model, bic, p, n, or NULL if no valid fit
#' @export
fit_DIAR_BIC <- function(y_h, y_base, F, k, sample_idx, max_p = 6, h = NULL, t_current = NULL) {
  best_bic <- Inf
  best_res <- NULL

  for (p in 0:max_p) {
    des <- build_design_DIAR(y_h, y_base, F, k, p, sample_idx, h = h, t_current = t_current)
    if (is.null(des)) next

    mod <- lm(des$y ~ des$X - 1)
    n <- length(des$y)
    k_param <- ncol(des$X)
    bic <- compute_bic(residuals(mod), n, k_param)

    if (bic < best_bic) {
      best_bic <- bic
      best_res <- list(model = mod, bic = bic, p = p, n = n)
    }
  }

  best_res
}

#' Predict DIAR at a given time origin
#'
#' @param fit_obj List returned from fit_DIAR_BIC
#' @param F Matrix: factor matrix
#' @param y_base Numeric vector: base series
#' @param t_idx Integer: time index for prediction
#' @param k Integer: number of factors to use
#' @return Numeric: predicted value, or NA if prediction fails
#' @export
predict_DIAR_at_origin <- function(fit_obj, F, y_base, t_idx, k) {
  if (is.null(fit_obj)) return(NA_real_)

  p <- fit_obj$p
  max_lag <- p

  # Need enough history for y-lags
  if (t_idx <= max_lag) return(NA_real_)

  # Factors at origin t_idx
  X_f <- F[t_idx, 1:k, drop = FALSE]

  # y-Lags
  if (p > 0) {
    y_lags <- sapply(1:p, function(j) y_base[t_idx - j])
    X_ylags <- matrix(y_lags, nrow = 1)
    colnames(X_ylags) <- paste0("ylag", 1:p)
    X_raw <- cbind(X_f, X_ylags)
  } else {
    X_raw <- X_f
  }

  if (any(!is.finite(c(X_raw)))) return(NA_real_)

  X_pred <- cbind(Intercept = 1, X_raw)
  beta   <- coef(fit_obj$model)

  # Safety: adjust length if lm did something unusual
  if (length(beta) != ncol(X_pred)) {
    beta <- beta[1:ncol(X_pred)]
  }

  as.numeric(X_pred %*% beta)
}

# ============================================================================
# DIAR-LAG (Direct forecast with factors and their lags, k <= 4)
# ============================================================================

#' Build DIAR-LAG design matrix
#'
#' @param y_h Numeric vector: h-step target series
#' @param y_base Numeric vector: base series for lags
#' @param F Matrix: factor matrix (T x Kmax)
#' @param k Integer: number of factors to use (must be <= 4 per Bae 2024)
#' @param p Integer: number of AR lags
#' @param m Integer: number of factor lags
#' @param sample_idx Integer vector: indices to use for estimation
#' @param h Integer: forecast horizon (for temporal constraint checking, optional)
#' @param t_current Integer: current time index (for temporal constraint checking, optional)
#' @return List with y and X, or NULL if insufficient data
build_design_DIARLAG <- function(y_h, y_base, F, k, p, m, sample_idx, h = NULL, t_current = NULL) {
  # TEMPORAL CONSTRAINT: Filter indices to respect t + h <= t_current
  if (!is.null(h) && !is.null(t_current)) {
    sample_idx <- sample_idx[sample_idx + h <= t_current]
  }

  max_lag <- max(p, m)
  idx <- sample_idx[sample_idx > max_lag & !is.na(y_h[sample_idx])]
  if (length(idx) == 0) return(NULL)

  y <- y_h[idx]

  # Factor lags (0 to m)
  X_f_list <- list()
  for (j in 0:m) {
    Fj <- F[idx - j, 1:k, drop = FALSE]
    colnames(Fj) <- paste0("F", 1:k, "_lag", j)
    X_f_list[[j + 1]] <- Fj
  }
  X_f <- do.call(cbind, X_f_list)

  # y-Lags
  if (p > 0) {
    X_ylags <- sapply(1:p, function(j) y_base[idx - j])
    colnames(X_ylags) <- paste0("ylag", 1:p)
    X_raw <- cbind(X_f, X_ylags)
  } else {
    X_raw <- X_f
  }

  ok <- stats::complete.cases(cbind(y, X_raw))
  y <- y[ok]
  X_raw <- X_raw[ok, , drop = FALSE]

  if (length(y) == 0) return(NULL)

  X <- cbind(Intercept = 1, X_raw)
  list(y = y, X = X)
}

#' Fit DIAR-LAG model with BIC-based lag selection
#'
#' @param y_h Numeric vector: h-step target series
#' @param y_base Numeric vector: base series for lags
#' @param F Matrix: factor matrix
#' @param k Integer: number of factors to use (must be <= 4)
#' @param sample_idx Integer vector: indices to use for estimation
#' @param max_p Maximum AR lag order to consider (default 6)
#' @param max_m Maximum factor lag order to consider (default 3)
#' @param h Integer: forecast horizon (for temporal constraint)
#' @param t_current Integer: current time index (for temporal constraint)
#' @return List with model, bic, p, m, n, or NULL if no valid fit
#' @export
fit_DIARLAG_BIC <- function(y_h, y_base, F, k, sample_idx,
                            max_p = 6, max_m = 3, h = NULL, t_current = NULL) {
  if (k > 4) {
    stop("DIAR-LAG allows at most k = 4 (per Bae 2024).")
  }

  best_bic <- Inf
  best_res <- NULL

  for (m in 1:max_m) {
    for (p in 0:max_p) {
      des <- build_design_DIARLAG(y_h, y_base, F, k, p, m, sample_idx, h = h, t_current = t_current)
      if (is.null(des)) next

      mod <- lm(des$y ~ des$X - 1)
      n <- length(des$y)
      k_param <- ncol(des$X)
      bic <- compute_bic(residuals(mod), n, k_param)

      if (bic < best_bic) {
        best_bic <- bic
        best_res <- list(model = mod, bic = bic, p = p, m = m, n = n)
      }
    }
  }

  best_res
}

#' Predict DIAR-LAG at a given time origin
#'
#' @param fit_obj List returned from fit_DIARLAG_BIC
#' @param F Matrix: factor matrix
#' @param y_base Numeric vector: base series
#' @param t_idx Integer: time index for prediction
#' @param k Integer: number of factors to use
#' @return Numeric: predicted value, or NA if prediction fails
#' @export
predict_DIARLAG_at_origin <- function(fit_obj, F, y_base, t_idx, k) {
  if (is.null(fit_obj)) return(NA_real_)

  p <- fit_obj$p
  m <- fit_obj$m
  max_lag <- max(p, m)

  # Need enough history for factor and y-lags
  if (t_idx <= max_lag) return(NA_real_)

  # Factor lags (0 to m)
  X_f_list <- list()
  for (j in 0:m) {
    Fj <- F[t_idx - j, 1:k, drop = FALSE]
    colnames(Fj) <- paste0("F", 1:k, "_lag", j)
    X_f_list[[j + 1]] <- Fj
  }
  X_f <- do.call(cbind, X_f_list)

  # y-Lags
  if (p > 0) {
    y_lags <- sapply(1:p, function(j) y_base[t_idx - j])
    X_ylags <- matrix(y_lags, nrow = 1)
    colnames(X_ylags) <- paste0("ylag", 1:p)
    X_raw <- cbind(X_f, X_ylags)
  } else {
    X_raw <- X_f
  }

  if (any(!is.finite(c(X_raw)))) return(NA_real_)

  X_pred <- cbind(Intercept = 1, X_raw)
  beta   <- coef(fit_obj$model)

  if (length(beta) != ncol(X_pred)) {
    beta <- beta[1:ncol(X_pred)]
  }

  as.numeric(X_pred %*% beta)
}
