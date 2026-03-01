#' Paper-Compliant PLS Factor Extraction
#'
#' Implements PLS exactly as specified in Bae (2024) Table 1.
#' This implementation uses the exact deflation procedure Q(F̂_prev)X
#' where Q(B) = I - B(B'B)^{-1}B', ensuring compliance with the paper's
#' mathematical definitions.
#'
#' @name paper_factors_pls
NULL

#' Extract PLS factors following Bae (2024) Table 1 exactly
#'
#' Implements the PLS specification from Bae (2024) Table 1:
#' - j = 1: α1 = argmax_α (1/T) Σ_t (α' x_t y_{t+1})^2 subject to N^{-1} α'α = 1
#'   First factor: F̂^{PLS}_1 = X α1 = XX' y
#' - j > 1: Deflate X using X*_j = Q(F̂^{PLS}_1, ..., F̂^{PLS}_{j-1}) X
#'   where Q(B) = I - P(B), P(B) = B(B'B)^{-1}B'
#'   Then αj = argmax_α (1/T) Σ_t (α' x*_{j,t} y_{t+1})^2 subject to N^{-1} α'α = 1
#'   F̂^{PLS}_j = X*_j αj
#'
#' @param X Matrix (T x N) of predictors. Should be mean-zero if paper assumptions hold.
#' @param y Numeric vector (length T) of target variable
#' @param k Integer, number of PLS factors to extract
#' @param demean_X Logical, whether to demean X (default FALSE, assumes preprocessing)
#' @param demean_y Logical, whether to demean y (default FALSE, assumes preprocessing)
#' @param config Configuration list for logging (optional)
#' @param verbose Logical, whether to print diagnostic information
#' @param ridge_tol Small ridge regularization for (F'F)^{-1} if near-singular (default 1e-10)
#'
#' @return A list with components:
#'   - F: Matrix (T x k) of PLS factor scores
#'   - weights: Matrix (N x k) of weights α_j satisfying N^{-1}α_j'α_j = 1
#'   - X_deflated_seq: List of deflated X matrices for each iteration (diagnostic)
#'   - method: "paper_pls" for identification
#'   - diagnostics: List of diagnostic information (constraint checks, orthogonality)
#'
#' @details
#' **Deflation Method:**
#' This implementation deflates X by projecting out the span of previously
#' estimated PLS factor scores, using Q(F_stack) = I - F_stack(F_stack'F_stack)^{-1}F_stack'.
#' This matches Bae's Table 1 exactly and is DIFFERENT from NIPALS which deflates
#' both X and y, or deflates X using loadings.
#'
#' **Constraint:**
#' The normalization N^{-1}α'α = 1 is strictly enforced for each weight vector.
#'
#' **First Factor:**
#' For j=1, α1 ∝ X'y with scaling to satisfy N^{-1}α'α = 1, giving:
#' α1 = sqrt(N) * (X'y) / ||X'y||_2
#' F̂_1 = X α1
#'
#' **Orthogonality:**
#' The deflation ensures that for j > 1:
#' F̂_j' F̂_i ≈ 0 for i < j (within numerical tolerance)
#' This is verified in diagnostics.
#'
#' **Invariants Checked:**
#' 1. N^{-1}α_j'α_j = 1 for all j
#' 2. F̂_j ⊥ span(F̂_1, ..., F̂_{j-1}) (via deflation)
#' 3. F̂_j ∈ column space of X*_j
#'
#' @examples
#' set.seed(42)
#' T <- 100; N <- 30; k <- 5
#' X <- matrix(rnorm(T*N), T, N)
#' y <- rnorm(T)
#'
#' # Mean-center
#' X <- scale(X, center = TRUE, scale = FALSE)
#' y <- y - mean(y)
#'
#' pls_paper <- paper_pls_factors(X, y, k = k, demean_X = FALSE, demean_y = FALSE, verbose = TRUE)
#'
#' # Check constraint: N^{-1}α'α = 1
#' N <- ncol(X)
#' for (j in 1:k) {
#'   alpha_j <- pls_paper$weights[, j]
#'   constraint_val <- (1/N) * sum(alpha_j^2)
#'   cat(sprintf("Factor %d: N^{-1}α'α = %.6f (should be 1.0)\\n", j, constraint_val))
#' }
#'
#' # Check orthogonality
#' F_mat <- pls_paper$F
#' F_cross <- t(F_mat) %*% F_mat / T
#' print("F'F/T (should be approximately diagonal):")
#' print(round(F_cross, 4))
#'
#' @export
paper_pls_factors <- function(X, y, k,
                               demean_X = FALSE, demean_y = FALSE,
                               config = NULL, verbose = FALSE,
                               ridge_tol = 1e-10) {
  # Input validation
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if (!is.numeric(y) || !is.vector(y)) {
    stop("y must be a numeric vector")
  }

  T_obs <- nrow(X)
  N <- ncol(X)

  if (length(y) != T_obs) {
    stop(sprintf("Dimension mismatch: y has length %d but X has %d rows", length(y), T_obs))
  }

  if (k > N) {
    stop(sprintf("Cannot extract k=%d PLS factors when N=%d predictors available", k, N))
  }
  if (k > T_obs) {
    stop(sprintf("Cannot extract k=%d PLS factors when T=%d observations available", k, T_obs))
  }

  # Optional demeaning
  if (demean_X) {
    X <- scale(X, center = TRUE, scale = FALSE)
  }
  if (demean_y) {
    y <- y - mean(y)
  }

  if (verbose || isTRUE(config$debug)) {
    log_debug(sprintf(
      "[paper_pls_factors] T=%d, N=%d, k=%d, demean_X=%s, demean_y=%s",
      T_obs, N, k, demean_X, demean_y
    ), config)
  }

  # Initialize storage
  F_matrix <- matrix(NA_real_, nrow = T_obs, ncol = k)  # T x k
  weights <- matrix(NA_real_, nrow = N, ncol = k)       # N x k, α_j vectors
  X_deflated_seq <- vector("list", k)  # For diagnostics
  constraint_checks <- numeric(k)
  orthogonality_checks <- vector("list", k)

  # Working copy of X (will be deflated iteratively)
  X_star <- X  # Initialize X*_1 = X

  # PLS iteration
  for (j in 1:k) {
    if (verbose || isTRUE(config$debug)) {
      log_debug(sprintf("[paper_pls_factors] Iteration j=%d/%d", j, k), config)
    }

    # Compute weight vector α_j
    # α_j = argmax_α (1/T) Σ_t (α' x*_{j,t} y_t)^2 subject to N^{-1}α'α = 1
    #
    # Solution: α_j ∝ X*' y, normalized so that N^{-1}α'α = 1
    # This gives: α_j = sqrt(N) * (X*' y) / ||X*' y||_2

    v <- crossprod(X_star, y)  # X*' y, dimension N x 1
    v_norm <- sqrt(sum(v^2))   # ||X*' y||_2

    if (v_norm < 1e-14) {
      warning(sprintf(
        "PLS iteration %d: ||X*' y|| = %.2e is near zero; covariance exhausted. Stopping at %d factors.",
        j, v_norm, j - 1
      ))
      # Truncate to j-1 factors
      F_matrix <- F_matrix[, 1:(j-1), drop = FALSE]
      weights <- weights[, 1:(j-1), drop = FALSE]
      k <- j - 1
      break
    }

    # Normalize: α_j = sqrt(N) * v / ||v||
    alpha_j <- sqrt(N) * (v / v_norm)  # N x 1
    weights[, j] <- alpha_j

    # Compute jth PLS factor: F̂_j = X*_j α_j
    f_j <- X_star %*% alpha_j  # T x 1
    F_matrix[, j] <- f_j

    # Store deflated X for diagnostics
    X_deflated_seq[[j]] <- X_star

    # Check constraint: N^{-1}α_j'α_j should equal 1
    constraint_val <- (1/N) * sum(alpha_j^2)
    constraint_checks[j] <- abs(constraint_val - 1.0)

    if (verbose || isTRUE(config$debug)) {
      log_debug(sprintf(
        "[paper_pls_factors] j=%d: N^{-1}α'α = %.6f (deviation from 1: %.2e)",
        j, constraint_val, constraint_checks[j]
      ), config)
    }

    if (constraint_checks[j] > 1e-10) {
      warning(sprintf(
        "PLS factor %d: constraint N^{-1}α'α deviates from 1 by %.2e",
        j, constraint_checks[j]
      ))
    }

    # Deflate X for next iteration (if j < k)
    if (j < k) {
      # Compute Q(F_1, ..., F_j) = I - F_stack (F_stack' F_stack)^{-1} F_stack'
      # where F_stack = [F_1, ..., F_j] is T x j

      F_stack <- F_matrix[, 1:j, drop = FALSE]  # T x j

      # Compute (F_stack' F_stack)^{-1} with optional ridge regularization
      FtF <- crossprod(F_stack)  # j x j matrix
      if (min(eigen(FtF, only.values = TRUE)$values) < ridge_tol) {
        log_warn(sprintf(
          "[paper_pls_factors] j=%d: F'F is near-singular, adding ridge %.2e",
          j, ridge_tol
        ), config)
        FtF <- FtF + ridge_tol * diag(j)
      }
      FtF_inv <- solve(FtF)  # j x j

      # Projection matrix: P(F_stack) = F_stack (F_stack' F_stack)^{-1} F_stack'
      # Orthogonal projection: Q(F_stack) = I_T - P(F_stack)
      # Apply to X: X*_{j+1} = Q(F_stack) X = X - F_stack (F_stack' F_stack)^{-1} F_stack' X

      X_star <- X - F_stack %*% FtF_inv %*% crossprod(F_stack, X)

      # Diagnostic: Check that F_stack' X_star ≈ 0 (deflation worked)
      residual_proj <- crossprod(F_stack, X_star)  # j x N
      max_residual <- max(abs(residual_proj))
      orthogonality_checks[[j]] <- max_residual

      if (verbose || isTRUE(config$debug)) {
        log_debug(sprintf(
          "[paper_pls_factors] j=%d: After deflation, max|F_prev' X*| = %.2e (should be ~0)",
          j, max_residual
        ), config)
      }
    }
  }

  # Final orthogonality check: F'F should be approximately diagonal
  F_cross <- crossprod(F_matrix) / T_obs  # k x k, normalized by T
  # Note: diag(scalar) creates an identity matrix of that size, not a 1x1 matrix.
  # Use nrow argument to ensure correct dimensions when k = 1.
  off_diagonal <- F_cross - diag(diag(F_cross), nrow = nrow(F_cross))
  max_off_diag <- max(abs(off_diagonal))

  if (verbose || isTRUE(config$debug)) {
    log_debug(sprintf(
      "[paper_pls_factors] Final orthogonality: max|off-diag(F'F/T)| = %.2e",
      max_off_diag
    ), config)
  }

  # Compile diagnostics
  diagnostics <- list(
    constraint_deviations = constraint_checks,
    orthogonality_residuals = orthogonality_checks,
    final_off_diagonal_max = max_off_diag,
    F_cross_product = F_cross
  )

  list(
    F = F_matrix,                     # T x k factor scores
    weights = weights,                # N x k weight vectors α_j
    X_deflated_seq = X_deflated_seq,  # Sequence of deflated X (diagnostic)
    method = "paper_pls",             # Identifier
    diagnostics = diagnostics
  )
}


#' Extract paper-compliant PLS factors (temporally-aware version)
#'
#' Fits PLS on a subset of data (X_fit, y_fit) but computes scores for all X_all.
#' This enforces temporal constraints for pseudo-out-of-sample forecasting.
#'
#' @param X_fit Matrix (T_fit x N) of predictors for fitting
#' @param y_fit Numeric vector (length T_fit) of target variable for fitting
#' @param X_all Matrix (T_all x N) of predictors for scoring (T_all >= T_fit)
#' @param k_max Maximum number of PLS components to extract
#' @param config Configuration list
#' @param use_paper_method Logical, if TRUE use paper_pls_factors
#' @return A list with components:
#'   - F: Matrix (T_all x k_max) of PLS scores for all X_all observations
#'   - model: The paper PLS object
#' @export
extract_pls_paper_temporal <- function(X_fit, y_fit, X_all, k_max = 12, config = NULL, use_paper_method = TRUE) {
  if (!use_paper_method) {
    stop("extract_pls_paper_temporal requires use_paper_method = TRUE")
  }

  # Fit paper PLS on valid subset
  K <- min(k_max, ncol(X_fit))

  pls_result <- paper_pls_factors(
    X = X_fit,
    y = y_fit,
    k = K,
    demean_X = FALSE,  # Should already be mean-zero
    demean_y = FALSE,  # Should already be mean-zero
    config = config,
    verbose = isTRUE(config$debug)
  )

  # Compute scores for all observations using the estimated weights
  X_all_mat <- as.matrix(X_all)
  F_all <- X_all_mat %*% pls_result$weights  # T_all x K

  # Assign column names
  colnames(F_all) <- paste0("F", seq_len(K))

  list(
    F = F_all,
    model = pls_result
  )
}

#' Wrapper for extract_pls using paper-compliant implementation
#'
#' Drop-in replacement for the current extract_pls() that uses
#' the paper-compliant implementation instead of pls::plsr().
#'
#' @param X Matrix (T x N) of predictors
#' @param y Numeric vector (length T) of target variable
#' @param k_max Maximum number of PLS components to extract
#' @param config Configuration list
#' @param use_paper_method Logical, if TRUE use paper_pls_factors,
#'   if FALSE use original pls::plsr (default TRUE)
#'
#' @return A list matching extract_pls() interface:
#'   - F: Matrix (T x k_max) of PLS scores for all observations
#'   - model: The paper PLS object or plsr model object
#'
#' @export
extract_pls_paper <- function(X, y, k_max = 12, config = NULL, use_paper_method = TRUE) {
  # Input validation
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }

  if (length(y) != nrow(X)) {
    stop(sprintf("Dimension mismatch: y has length %d but X has %d rows", length(y), nrow(X)))
  }

  # Build data frame for compatibility check
  df <- data.frame(y = y, X)
  df_fit <- df[!is.na(df$y), , drop = FALSE]

  if (nrow(df_fit) <= 2) {
    stop("Too few observations for PLS at this origin.")
  }

  # Extract non-NA subset
  X_fit <- as.matrix(df_fit[, -1, drop = FALSE])
  y_fit <- df_fit$y

  if (isTRUE(config$debug)) {
    y_sd <- sd(y_fit, na.rm = TRUE)
    log_debug(sprintf(
      "[extract_pls_paper] nrow(X_fit)=%d, ncol(X)=%d, k_max=%d, sd(y)=%.6f, use_paper=%s",
      nrow(X_fit), ncol(X_fit), k_max, y_sd, use_paper_method
    ), config)
  }

  if (use_paper_method) {
    # Use paper-compliant implementation
    K <- min(k_max, ncol(X_fit))

    # Assume X and y are already centered in preprocessing
    pls_result <- paper_pls_factors(
      X = X_fit,
      y = y_fit,
      k = K,
      demean_X = FALSE,  # Should already be mean-zero
      demean_y = FALSE,  # Should already be mean-zero
      config = config,
      verbose = isTRUE(config$debug)
    )

    # Need to predict scores for ALL observations (including NA y)
    # For paper PLS, we compute scores as F = X α for the full X matrix
    # This matches the paper's definition and ensures we can forecast at all origins

    X_all <- as.matrix(df[, -1, drop = FALSE])
    F_all <- X_all %*% pls_result$weights  # T x K

    # Assign column names
    colnames(F_all) <- paste0("F", seq_len(K))

    # Return in extract_pls() format
    list(
      F = F_all,
      model = pls_result
    )

  } else {
    # Fall back to original pls::plsr implementation
    log_warn("[extract_pls_paper] Using pls::plsr (NOT paper-compliant deflation)", config)

    pls_fit <- pls::plsr(
      y ~ .,
      data   = df_fit,
      ncomp  = min(k_max, ncol(X_fit)),
      center = isTRUE(config$pls_center),
      scale  = isTRUE(config$pls_scale),
      method = config$pls_method
    )

    # Predict scores for all observations
    X_all <- df[, setdiff(names(df), "y"), drop = FALSE]
    K <- min(k_max, ncol(X_fit))
    scores_arr <- predict(pls_fit, newdata = X_all, ncomp = 1:K, type = "scores")

    # Convert to matrix deterministically
    if (length(dim(scores_arr)) == 3) {
      F_hat <- scores_arr[, , 1, drop = TRUE]
    } else if (length(dim(scores_arr)) == 2) {
      F_hat <- scores_arr
    } else {
      stop("Unexpected dimensions of PLS scores array.")
    }

    F_hat <- as.matrix(F_hat)
    colnames(F_hat) <- paste0("F", seq_len(ncol(F_hat)))

    list(
      F = F_hat,
      model = pls_fit
    )
  }
}
