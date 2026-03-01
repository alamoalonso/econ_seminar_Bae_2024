#' Factor Number Selection (k) Decision Rules
#'
#' Pluggable decision rules for selecting the number of factors k.
#' Implements Bai-Ng BIC3 for PCA and Onatski (2010) for PLS.
#'
#' @name k_selection
NULL

# ============================================================================
# BN-BIC (Bai-Ng BIC3) for PCA
# ============================================================================

#' Compute BN-BIC criterion and select k for PCA
#'
#' Implements BIC3 from Bai & Ng (2002) for selecting the number of PCA factors.
#' At each origin, computes reconstruction residual variance V(k) for k=1..k_max
#' and selects k_hat = argmin BIC3(k).
#'
#' @param X Matrix (T x N) of predictors (should be centered/scaled as per repo convention)
#' @param k_max Maximum k to consider
#' @param sigma_sq_rule How to compute sigma^2: "v_kmax" (default)
#' @param min_k Minimum k allowed (default: 1)
#' @param config Configuration list for logging
#' @return List with:
#'   - k_hat: Selected number of factors (>= min_k)
#'   - BIC_values: Vector of BIC3 values for k = 1..k_max
#'   - V_values: Vector of V(k) for k = 1..k_max
#'   - eigenvalues: Eigenvalues from decomposition
#'   - warning: Warning message if k_hat was adjusted to min_k
#' @export
#'
#' @details
#' Formula (BIC3 from Bai & Ng 2002):
#' BIC3(k) = ln(V(k)) + k * sigma^2 * ((N + T - k) / (N*T)) * ln(N*T)
#'
#' where V(k) = (1/(N*T)) * sum_{i,t} (x_it - lambda_i' F_t)^2 is the
#' mean squared reconstruction error using k factors.
#'
#' sigma^2 = V(k_max) by default (common convention).
compute_bn_bic_k <- function(X, k_max, sigma_sq_rule = "v_kmax", min_k = 1L, config = NULL) {

  if (!is.matrix(X)) X <- as.matrix(X)

  T_obs <- nrow(X)
  N <- ncol(X)

  # Validate inputs
  if (k_max < 1) stop("k_max must be >= 1")
  if (k_max > min(T_obs, N)) {
    k_max <- min(T_obs, N)
    log_warn(sprintf("[BN-BIC] k_max reduced to min(T,N) = %d", k_max), config)
  }

  # Step 1: Compute SVD once (efficient)
  # X = U * D * V'
  # PCA: F = U * D (scores), Lambda = V (loadings)
  # Reconstruction with k factors: X_k = U[,1:k] * D[1:k] * V[,1:k]'
  # Residual: E_k = X - X_k

  # Center X if not already centered (BN-BIC assumes centered data)
  X_centered <- scale(X, center = TRUE, scale = FALSE)

  svd_result <- svd(X_centered, nu = k_max, nv = k_max)
  d <- svd_result$d
  U <- svd_result$u
  V <- svd_result$v

  # Total sum of squares (for k=0 reference)
  total_ss <- sum(X_centered^2)

  # Step 2: Compute V(k) for k = 1..k_max using cumulative explained variance
  # V(k) = (1/(N*T)) * ||X - X_k||^2 = (1/(N*T)) * (total_ss - sum(d[1:k]^2))
  # This is efficient: no need to reconstruct X_k explicitly

  cumsum_d2 <- cumsum(d[1:k_max]^2)
  residual_ss <- total_ss - cumsum_d2  # Residual SS for k = 1..k_max
  V_values <- residual_ss / (N * T_obs)

  # Handle numerical issues (V(k) should be non-negative)
  V_values <- pmax(V_values, .Machine$double.eps)

  # Step 3: Compute sigma^2
  if (sigma_sq_rule == "v_kmax") {
    sigma_sq <- V_values[k_max]
  } else {
    # Default fallback
    sigma_sq <- V_values[k_max]
  }

  # Step 4: Compute BIC3 for k = 1..k_max
  # BIC3(k) = ln(V(k)) + k * sigma^2 * ((N + T - k) / (N*T)) * ln(N*T)
  k_vec <- 1:k_max
  penalty <- k_vec * sigma_sq * ((N + T_obs - k_vec) / (N * T_obs)) * log(N * T_obs)
  BIC_values <- log(V_values) + penalty

  # Step 5: Select k_hat = argmin BIC3(k)
  k_hat_raw <- which.min(BIC_values)

  # Step 6: Enforce k_hat >= min_k
  warning_msg <- NULL
  if (k_hat_raw < min_k) {
    warning_msg <- sprintf("BN-BIC selected k=%d but min_k=%d; setting k_hat=%d",
                           k_hat_raw, min_k, min_k)
    k_hat <- min_k
  } else {
    k_hat <- k_hat_raw
  }

  list(
    k_hat = as.integer(k_hat),
    k_hat_raw = as.integer(k_hat_raw),
    BIC_values = BIC_values,
    V_values = V_values,
    eigenvalues = d^2,  # Eigenvalues of X'X (proportional to variance explained)
    sigma_sq = sigma_sq,
    T_obs = T_obs,
    N = N,
    warning = warning_msg
  )
}


# ============================================================================
# Onatski (2010) for PLS
# ============================================================================

#' Compute eigenvalues of sample covariance matrix
#'
#' Computes eigenvalues of (1/T) * X' * X (N x N covariance in variable dimension).
#' Used by Onatski (2010) edge detection method.
#'
#' @param X Matrix (T x N) of predictors (should be preprocessed as per PLS convention)
#' @return Vector of eigenvalues in descending order
#' @export
compute_cov_eigenvalues <- function(X) {
  if (!is.matrix(X)) X <- as.matrix(X)

  T_obs <- nrow(X)

  # Covariance matrix in variable dimension: (1/T) * X' * X
  # This is N x N
  cov_matrix <- crossprod(X) / T_obs

  # Eigenvalues only (we don't need eigenvectors for Onatski)
  eigenvalues <- eigen(cov_matrix, symmetric = TRUE, only.values = TRUE)$values

  # Sort descending (should already be, but ensure)
  sort(eigenvalues, decreasing = TRUE)
}


#' Compute Onatski (2010) criterion and select k
#'
#' Implements the Onatski (2010) eigenvalue edge detection method for
#' selecting the number of factors. Uses eigenvalues of the sample covariance
#' of X (unsupervised, not dependent on y).
#'
#' @param X Matrix (T x N) of predictors (preprocessed consistently with PLS)
#' @param r_max Maximum r for edge detection (default: 12)
#' @param min_k Minimum k allowed (default: 1)
#' @param k_max Maximum k allowed (default: NULL, no upper bound)
#' @param config Configuration list for logging
#' @return List with:
#'   - k_hat: Selected number of factors (min_k <= k_hat <= k_max)
#'   - eigenvalues: All eigenvalues (descending order)
#'   - u_hat: Edge detection threshold
#'   - delta: Shrinking threshold
#'   - r_max_used: Actual r_max used (may be adjusted)
#'   - warning: Warning message if k_hat was adjusted
#' @export
#'
#' @details
#' Onatski (2010) method:
#' 1. Compute eigenvalues lambda_1 >= ... >= lambda_N of sample covariance
#' 2. Determine r_max such that 2*r_max + 1 <= min(N, T)
#' 3. Compute w = 2^(2/3) / (2^(2/3) - 1) approx 2.7321
#' 4. Compute u_hat = w * lambda_{r_max+1} + (1-w) * lambda_{2*r_max+1}
#' 5. Compute delta = max(N^{-2/5}, T^{-2/5})
#' 6. Select k_hat = #{i : lambda_i > (1 + delta) * u_hat}, clamped to [min_k, k_max]
compute_onatski_k <- function(X, r_max = 12L, min_k = 1L, k_max = NULL, config = NULL) {

  if (!is.matrix(X)) X <- as.matrix(X)

  T_obs <- nrow(X)
  N <- ncol(X)

  # Step 1: Compute eigenvalues
  eigenvalues <- compute_cov_eigenvalues(X)

  # Step 2: Adjust r_max if constraint violated
  # Constraint: 2*r_max + 1 <= min(N, T)
  # => r_max <= (min(N, T) - 1) / 2
  r_max_constraint <- floor((min(N, T_obs) - 1) / 2)
  r_max_used <- min(r_max, r_max_constraint)

  warning_msg <- NULL

  if (r_max_used < 1) {
    # Cannot apply Onatski: too few observations/variables
    warning_msg <- sprintf(
      "Onatski infeasible: min(N=%d, T=%d) too small for r_max >= 1; returning k_hat=%d",
      N, T_obs, min_k
    )
    log_warn(sprintf("[Onatski] %s", warning_msg), config)

    return(list(
      k_hat = as.integer(min_k),
      k_hat_raw = 0L,
      eigenvalues = eigenvalues,
      u_hat = NA_real_,
      delta = NA_real_,
      r_max_used = 0L,
      warning = warning_msg
    ))
  }

  if (r_max_used < r_max) {
    log_warn(sprintf("[Onatski] r_max reduced from %d to %d due to constraint 2*r_max+1 <= min(N,T)",
                     r_max, r_max_used), config)
  }

  # Step 3: Compute w
  w <- 2^(2/3) / (2^(2/3) - 1)  # approx 2.7321

  # Step 4: Compute u_hat
  # Need eigenvalues at indices r_max+1 and 2*r_max+1
  idx1 <- r_max_used + 1
  idx2 <- 2 * r_max_used + 1

  if (idx2 > length(eigenvalues)) {
    # This shouldn't happen given the constraint, but safety check
    warning_msg <- sprintf(
      "Onatski: eigenvalue index %d exceeds available %d; returning k_hat=%d",
      idx2, length(eigenvalues), min_k
    )
    log_warn(sprintf("[Onatski] %s", warning_msg), config)

    return(list(
      k_hat = as.integer(min_k),
      k_hat_raw = 0L,
      eigenvalues = eigenvalues,
      u_hat = NA_real_,
      delta = NA_real_,
      r_max_used = r_max_used,
      warning = warning_msg
    ))
  }

  lambda_r1 <- eigenvalues[idx1]
  lambda_r2 <- eigenvalues[idx2]

  u_hat <- w * lambda_r1 + (1 - w) * lambda_r2

  # Step 5: Compute delta
  delta <- max(N^(-2/5), T_obs^(-2/5))

  # Step 6: Select k_hat = #{i : lambda_i > (1 + delta) * u_hat}
  threshold <- (1 + delta) * u_hat
  k_hat_raw <- sum(eigenvalues > threshold)

  # Step 7: Enforce min_k <= k_hat <= k_max
  k_hat <- k_hat_raw

  if (k_hat < min_k) {
    if (is.null(warning_msg)) {
      warning_msg <- sprintf("Onatski selected k=%d but min_k=%d; setting k_hat=%d",
                             k_hat_raw, min_k, min_k)
    }
    k_hat <- min_k
  }

  if (!is.null(k_max) && k_hat > k_max) {
    if (is.null(warning_msg)) {
      warning_msg <- sprintf("Onatski selected k=%d but k_max=%d; setting k_hat=%d",
                             k_hat_raw, k_max, k_max)
    } else {
      warning_msg <- paste(warning_msg, sprintf("; also clamped to k_max=%d", k_max))
    }
    k_hat <- k_max
  }

  list(
    k_hat = as.integer(k_hat),
    k_hat_raw = as.integer(k_hat_raw),
    eigenvalues = eigenvalues,
    u_hat = u_hat,
    delta = delta,
    threshold = threshold,
    r_max_used = r_max_used,
    T_obs = T_obs,
    N = N,
    warning = warning_msg
  )
}


# ============================================================================
# Unified k-Selection Dispatcher
# ============================================================================

#' Select k using specified decision rule
#'
#' Dispatcher function that calls the appropriate k-selection method.
#'
#' @param X Matrix (T x N) of predictors
#' @param k_rule Rule name: "bn_bic" or "onatski"
#' @param k_max Upper bound for k
#' @param config Configuration list (includes k_selection_settings)
#' @param context List with contextual info for logging (series_id, origin_index, origin_date)
#' @return List with k_hat and rule-specific details
#' @export
select_k_dynamic <- function(X, k_rule, k_max, config = NULL, context = NULL) {

  # Get settings
  settings <- config$k_selection_settings
  min_k <- if (!is.null(settings$min_k)) settings$min_k else 1L

  # Dispatch to appropriate method
  if (k_rule == "bn_bic") {
    sigma_sq_rule <- if (!is.null(settings$bn_bic_sigma_sq)) {
      settings$bn_bic_sigma_sq
    } else {
      "v_kmax"
    }
    result <- compute_bn_bic_k(X, k_max, sigma_sq_rule, min_k, config)

  } else if (k_rule == "onatski") {
    r_max <- if (!is.null(settings$onatski_r_max)) {
      settings$onatski_r_max
    } else {
      12L
    }
    result <- compute_onatski_k(X, r_max, min_k, k_max, config)

  } else {
    stop(sprintf("Unknown k_rule: '%s'. Must be 'bn_bic' or 'onatski'.", k_rule))
  }

  # Log warning if k was adjusted to min_k
  if (!is.null(result$warning) && !is.null(context)) {
    log_warn(sprintf(
      "[k-selection] %s | series=%s, origin_idx=%s, origin_date=%s",
      result$warning,
      context$series_id %||% "?",
      context$origin_index %||% "?",
      as.character(context$origin_date %||% "?")
    ), config)
  } else if (!is.null(result$warning)) {
    log_warn(sprintf("[k-selection] %s", result$warning), config)
  }

  # Add rule info
  result$k_rule <- k_rule

  result
}

# Helper: null-coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x
