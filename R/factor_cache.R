#' Factor Extraction Cache
#'
#' Per-origin caching for factor extraction to enable efficient evaluation
#' of both grid and dynamic k-selection methods without redundant computation.
#'
#' @name factor_cache
NULL

#' Extract and cache factors at a single origin for all specs
#'
#' Extracts PCA and/or PLS factors once per origin, caching the results
#' for efficient use by both grid and dynamic specs.
#'
#' @param X_train Matrix (T_train x N) of training predictors (up to origin)
#' @param y_train Numeric vector of training target (for PLS; may be shorter due to h offset)
#' @param X_full Matrix of all predictors up to origin (for scoring)
#' @param valid_pls_idx Indices in X_full that are valid for PLS fitting
#' @param factor_specs List of factor_spec objects
#' @param config Configuration list
#' @return List with cached extraction results:
#'   - pca: List with F, loadings, eigenvalues, V_values (if any PCA spec)
#'   - pls: List with F, model, X_pls (preprocessed X for Onatski)
#' @export
extract_factors_for_specs <- function(X_train, y_train, X_full, valid_pls_idx,
                                       factor_specs, config = NULL) {
  cache <- list(pca = NULL, pls = NULL)

  # Determine max k needed for each factor method
  k_max_pca <- get_max_k_for_method(factor_specs, "PCA")
  k_max_pls <- max(
    get_max_k_for_method(factor_specs, "PLS"),
    get_max_k_for_method(factor_specs, "1-PLS")
  )

  # Check if we need BN-BIC (requires V_values from PCA)
  needs_bn_bic <- has_dynamic_rule(factor_specs, "bn_bic")

  # Check if we need Onatski (requires eigenvalues from PLS-preprocessed X)
  needs_onatski <- has_dynamic_rule(factor_specs, "onatski")

  # ---- PCA Extraction ----
  if (k_max_pca > 0) {
    cache$pca <- extract_pca_cached(
      X_train = X_train,
      k_max = k_max_pca,
      compute_V_k = needs_bn_bic,
      config = config
    )
  }

  # ---- PLS Extraction ----
  if (k_max_pls > 0 && length(valid_pls_idx) >= 2) {
    # Get PLS preprocessing settings
    pls_center <- isTRUE(config$pls_center)
    pls_scale <- isTRUE(config$pls_scale)

    cache$pls <- extract_pls_cached(
      X_full = X_full,
      y_train = y_train,
      valid_pls_idx = valid_pls_idx,
      k_max = k_max_pls,
      pls_center = pls_center,
      pls_scale = pls_scale,
      compute_eigenvalues = needs_onatski,
      config = config
    )
  }

  cache
}


#' Extract PCA factors with caching support
#'
#' Extracts PCA factors and optionally computes V(k) for BN-BIC.
#' Uses the same PCA implementation as the repo's standard extraction.
#'
#' @param X_train Matrix (T x N) of training predictors
#' @param k_max Maximum number of components
#' @param compute_V_k Logical: compute V(k) values for BN-BIC?
#' @param config Configuration list
#' @return List with F, loadings, eigenvalues, V_values (if computed)
#' @keywords internal
extract_pca_cached <- function(X_train, k_max, compute_V_k = FALSE, config = NULL) {

  if (!is.matrix(X_train)) X_train <- as.matrix(X_train)

  T_obs <- nrow(X_train)
  N <- ncol(X_train)

  # Adjust k_max if needed
  k_max <- min(k_max, T_obs, N)

  # Use config settings for centering/scaling
  pca_center <- isTRUE(config$pca_center)
  pca_scale <- isTRUE(config$pca_scale)

  # Dispatch to paper-compliant or library implementation
  if (isTRUE(config$use_paper_pca)) {
    # Paper PCA uses custom normalization
    pca_out <- extract_pca_paper(X_train, k_max = k_max, config = config, use_paper_method = TRUE)
  } else {
    # Library PCA (prcomp)
    pca_out <- extract_pca(X_train, k_max = k_max, config = config)
  }

  result <- list(
    F = pca_out$F,
    loadings = pca_out$loadings,
    eigenvalues = NULL,
    V_values = NULL
  )

  # Compute eigenvalues and V(k) if needed for BN-BIC
  if (compute_V_k) {
    # Get eigenvalues from SVD of centered X
    X_centered <- scale(X_train, center = pca_center, scale = pca_scale)
    svd_result <- svd(X_centered, nu = 0, nv = 0)  # Only need singular values
    d <- svd_result$d

    result$eigenvalues <- d^2

    # Compute V(k) = (1/(N*T)) * ||X - X_k||^2
    total_ss <- sum(X_centered^2)
    k_max_svd <- min(k_max, length(d))
    cumsum_d2 <- cumsum(d[1:k_max_svd]^2)
    residual_ss <- total_ss - cumsum_d2
    result$V_values <- pmax(residual_ss / (N * T_obs), .Machine$double.eps)
  }

  result
}


#' Extract PLS factors with caching support
#'
#' Extracts PLS factors with temporal constraint handling and optionally
#' computes eigenvalues for Onatski from the same preprocessed X.
#'
#' @param X_full Matrix of all predictors up to origin
#' @param y_train Numeric vector of training target
#' @param valid_pls_idx Indices valid for PLS fitting (respecting t+h <= t_origin)
#' @param k_max Maximum number of components
#' @param pls_center Logical: center predictors?
#' @param pls_scale Logical: scale predictors?
#' @param compute_eigenvalues Logical: compute eigenvalues for Onatski?
#' @param config Configuration list
#' @return List with F, model, eigenvalues (if computed), X_pls_preprocessed
#' @keywords internal
extract_pls_cached <- function(X_full, y_train, valid_pls_idx, k_max,
                                pls_center = TRUE, pls_scale = FALSE,
                                compute_eigenvalues = FALSE, config = NULL) {

  if (!is.matrix(X_full)) X_full <- as.matrix(X_full)

  # Extract valid training data
  X_pls_fit <- X_full[valid_pls_idx, , drop = FALSE]
  y_pls_fit <- y_train[valid_pls_idx]

  # Remove NA y values
  valid_y <- !is.na(y_pls_fit)
  X_pls_fit <- X_pls_fit[valid_y, , drop = FALSE]
  y_pls_fit <- y_pls_fit[valid_y]

  result <- list(
    F = NULL,
    model = NULL,
    eigenvalues = NULL,
    X_pls_preprocessed = NULL
  )

  if (nrow(X_pls_fit) < 2) {
    log_warn("[PLS Cache] Too few observations for PLS extraction", config)
    result$F <- matrix(NA_real_, nrow = nrow(X_full), ncol = k_max)
    return(result)
  }

  k_max <- min(k_max, nrow(X_pls_fit), ncol(X_pls_fit))

  # Dispatch to paper-compliant or library implementation
  if (isTRUE(config$use_paper_pls)) {
    pls_out <- extract_pls_paper_temporal(
      X_fit = X_pls_fit,
      y_fit = y_pls_fit,
      X_all = X_full,
      k_max = k_max,
      config = config,
      use_paper_method = TRUE
    )
  } else {
    pls_out <- extract_pls_temporal(
      X_fit = X_pls_fit,
      y_fit = y_pls_fit,
      X_all = X_full,
      k_max = k_max,
      config = config
    )
  }

  result$F <- pls_out$F
  result$model <- pls_out$model

  # Compute eigenvalues for Onatski from PLS-preprocessed X
  # IMPORTANT: Use same preprocessing as PLS to ensure consistency
  if (compute_eigenvalues) {
    # Apply same center/scale as PLS
    X_preprocessed <- scale(X_pls_fit, center = pls_center, scale = pls_scale)
    result$eigenvalues <- compute_cov_eigenvalues(X_preprocessed)
    result$X_pls_preprocessed <- X_preprocessed
  }

  result
}


#' Generate forecasts for a single spec at one origin
#'
#' Uses cached factors to generate forecasts for a factor_spec.
#' For grid specs: generates forecasts for k=1..k_max
#' For dynamic specs: computes k_hat and generates single forecast
#'
#' @param spec Factor specification
#' @param cache Cached factor extraction results
#' @param y_h Target variable (h-step ahead)
#' @param y_base Base target for AR lags
#' @param idx_rec Indices for recursive scheme
#' @param idx_roll Indices for rolling scheme
#' @param t_idx Current time index
#' @param h Horizon
#' @param config Configuration list
#' @param context Contextual info for logging
#' @return List with forecasts for this spec
#' @keywords internal
generate_forecasts_for_spec <- function(spec, cache, y_h, y_base,
                                         idx_rec, idx_roll, t_idx, h,
                                         config = NULL, context = NULL) {

  factor_method <- spec$factor_method
  k_mode <- spec$k_mode
  k_rule <- spec$k_rule
  k_max_spec <- spec$k_max

  # Select appropriate factor matrix
  if (factor_method == "PCA") {
    if (is.null(cache$pca)) {
      return(list(error = "PCA cache not available"))
    }
    F_used <- cache$pca$F
    eigenvalues <- cache$pca$eigenvalues
    V_values <- cache$pca$V_values
  } else {
    # PLS or 1-PLS
    if (is.null(cache$pls)) {
      return(list(error = "PLS cache not available"))
    }
    F_used <- cache$pls$F
    eigenvalues <- cache$pls$eigenvalues
  }

  # Ensure k_max doesn't exceed available factors
  k_max_used <- min(k_max_spec, ncol(F_used))

  result <- list(
    spec_id = spec$id,
    k_mode = k_mode,
    k_rule = k_rule,
    forecasts = list()
  )

  if (k_mode == "grid") {
    # Grid: generate forecasts for k = 1..k_max
    result$k_values <- 1:k_max_used
    result$k_hat <- NA_integer_  # Not applicable for grid

    for (k in 1:k_max_used) {
      result$forecasts[[paste0("k", k)]] <- list(
        k = k,
        DI_rec = NA_real_,
        DI_roll = NA_real_,
        DIAR_rec = NA_real_,
        DIAR_roll = NA_real_,
        DLAG_rec = NA_real_,
        DLAG_roll = NA_real_
      )
    }

  } else {
    # Dynamic: compute k_hat and generate single forecast
    if (factor_method == "PCA" && k_rule == "bn_bic") {
      # Use cached eigenvalues/V_values if available, otherwise recompute
      if (!is.null(V_values)) {
        # Recompute BIC from cached V_values
        T_obs <- nrow(F_used)
        N <- ncol(cache$pca$loadings)
        sigma_sq <- V_values[k_max_used]

        k_vec <- 1:k_max_used
        penalty <- k_vec * sigma_sq * ((N + T_obs - k_vec) / (N * T_obs)) * log(N * T_obs)
        BIC_values <- log(V_values[1:k_max_used]) + penalty

        k_hat_raw <- which.min(BIC_values)
        min_k <- config$k_selection_settings$min_k %||% 1L
        k_hat <- max(k_hat_raw, min_k)

        if (k_hat_raw < min_k && !is.null(context)) {
          log_warn(sprintf(
            "[BN-BIC] k_hat_raw=%d < min_k=%d; using k_hat=%d | series=%s, origin=%d",
            k_hat_raw, min_k, k_hat, context$series_id, t_idx
          ), config)
        }
      } else {
        # Fallback: should not happen if cache is properly computed
        k_hat <- k_max_used
        log_warn("[BN-BIC] V_values not cached; using k_max as fallback", config)
      }

    } else if (k_rule == "onatski") {
      # Use cached eigenvalues from PLS-preprocessed X
      if (!is.null(eigenvalues)) {
        settings <- config$k_selection_settings
        r_max <- settings$onatski_r_max %||% 12L
        min_k <- settings$min_k %||% 1L

        # Onatski computation from eigenvalues
        on_result <- compute_onatski_from_eigenvalues(eigenvalues, r_max, min_k, config)
        k_hat <- on_result$k_hat

        if (on_result$k_hat_raw < min_k && !is.null(context)) {
          log_warn(sprintf(
            "[Onatski] k_hat_raw=%d < min_k=%d; using k_hat=%d | series=%s, origin=%d",
            on_result$k_hat_raw, min_k, k_hat, context$series_id, t_idx
          ), config)
        }
      } else {
        # Fallback
        k_hat <- k_max_used
        log_warn("[Onatski] eigenvalues not cached; using k_max as fallback", config)
      }
    } else {
      stop(sprintf("Unknown k_rule '%s' for factor_method '%s'", k_rule, factor_method))
    }

    # Ensure k_hat <= k_max_used
    k_hat <- min(k_hat, k_max_used)

    result$k_values <- k_hat
    result$k_hat <- k_hat
    result$forecasts[[paste0("k", k_hat)]] <- list(
      k = k_hat,
      DI_rec = NA_real_,
      DI_roll = NA_real_,
      DIAR_rec = NA_real_,
      DIAR_roll = NA_real_,
      DLAG_rec = NA_real_,
      DLAG_roll = NA_real_
    )
  }

  result
}


#' Compute Onatski k from pre-computed eigenvalues
#'
#' @param eigenvalues Vector of eigenvalues in descending order
#' @param r_max Maximum r for edge detection
#' @param min_k Minimum k allowed
#' @param config Configuration list
#' @return List with k_hat, k_hat_raw, etc.
#' @keywords internal
compute_onatski_from_eigenvalues <- function(eigenvalues, r_max = 12L, min_k = 1L, config = NULL) {

  N <- length(eigenvalues)

  # We don't have T_obs directly from eigenvalues, but for Onatski's delta

  # we use the original N (from covariance matrix dimension)
  # The constraint check was done during extraction

  # Adjust r_max
  r_max_constraint <- floor((N - 1) / 2)
  r_max_used <- min(r_max, r_max_constraint)

  if (r_max_used < 1) {
    return(list(
      k_hat = as.integer(min_k),
      k_hat_raw = 0L,
      warning = "r_max_used < 1"
    ))
  }

  # Compute w
  w <- 2^(2/3) / (2^(2/3) - 1)

  # Compute u_hat
  idx1 <- r_max_used + 1
  idx2 <- 2 * r_max_used + 1

  if (idx2 > N) {
    return(list(
      k_hat = as.integer(min_k),
      k_hat_raw = 0L,
      warning = "eigenvalue index out of bounds"
    ))
  }

  u_hat <- w * eigenvalues[idx1] + (1 - w) * eigenvalues[idx2]

  # Use N for delta (this is an approximation; full version would need T_obs)
  delta <- N^(-2/5)

  # Select k
  threshold <- (1 + delta) * u_hat
  k_hat_raw <- sum(eigenvalues > threshold)

  k_hat <- max(k_hat_raw, min_k)

  list(
    k_hat = as.integer(k_hat),
    k_hat_raw = as.integer(k_hat_raw),
    u_hat = u_hat,
    delta = delta,
    threshold = threshold
  )
}

# Helper: null-coalescing operator (if not already defined)
if (!exists("%||%")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}
