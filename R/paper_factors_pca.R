#' Paper-Compliant PCA Factor Extraction
#'
#' Implements PCA exactly as specified in Bae (2024) Section 2.2.
#' This implementation ensures deterministic compliance with the paper's
#' mathematical definitions, including the specific normalization N^{-1}Λ'Λ = I_r
#' and factor formula F̂ = N^{-1} X Λ̂.
#'
#' @name paper_factors_pca
NULL

#' Extract PCA factors following Bae (2024) exactly
#'
#' Implements the PCA specification from Bae (2024) Section 2.2:
#' - minimize (1/(NT)) Σ_t (x_t − Λ F_t)'(x_t − Λ F_t)
#' - subject to N^{-1} Λ'Λ = I_r
#' - Solution: Λ̂_PCA = scaled r eigenvectors of Σ̂_X = (1/T) X'X
#' - Factors: F̂_PCA = N^{-1} X Λ̂_PCA
#'
#' @param X Matrix (T x N) of predictors. Should already be demeaned if the
#'   paper's assumption of mean-zero data is satisfied. If X is not centered,
#'   set demean=TRUE.
#' @param r Integer, number of factors to extract (k_max in forecasting context)
#' @param demean Logical, whether to demean X before computing covariance.
#'   Default FALSE assumes X is already mean-zero per paper preprocessing.
#' @param config Configuration list for logging (optional)
#' @param verbose Logical, whether to print diagnostic information
#'
#' @return A list with components:
#'   - F: Matrix (T x r) of factor scores, F̂ = N^{-1} X Λ̂
#'   - Lambda: Matrix (N x r) of factor loadings satisfying N^{-1}Λ'Λ = I_r
#'   - eigenvalues: Vector of top r eigenvalues of Σ̂_X
#'   - Sigma_X: Matrix (N x N) sample covariance Σ̂_X = (1/T) X'X
#'   - method: "paper_pca" for identification
#'
#' @details
#' **Normalization:**
#' Unlike base R's prcomp() which returns Λ with Λ'Λ = I_r,
#' this function enforces Bae's normalization: N^{-1} Λ'Λ = I_r.
#' This means columns of Λ have norm √N, not 1.
#'
#' **Factor Computation:**
#' Factors are computed as F̂ = N^{-1} X Λ̂, not simply X Λ̂.
#' This N^{-1} scaling is critical for matching the paper's specification.
#'
#' **Invariants Verified:**
#' The function checks that N^{-1} Λ'Λ ≈ I_r within tolerance.
#' If verbose=TRUE, reports deviation from identity.
#'
#' @examples
#' set.seed(123)
#' T <- 100; N <- 50; r <- 5
#' X <- matrix(rnorm(T*N), T, N)
#' X <- scale(X, center = TRUE, scale = FALSE)  # Mean-zero
#'
#' pca_paper <- paper_pca_factors(X, r = r, demean = FALSE, verbose = TRUE)
#'
#' # Verify normalization
#' N <- ncol(X)
#' Lambda <- pca_paper$Lambda
#' I_r <- diag(r)
#' norm_check <- (1/N) * t(Lambda) %*% Lambda
#' print(max(abs(norm_check - I_r)))  # Should be near zero
#'
#' # Verify factor formula
#' F_check <- (1/N) * X %*% Lambda
#' print(max(abs(F_check - pca_paper$F)))  # Should be near zero
#'
#' @export
paper_pca_factors <- function(X, r, demean = FALSE, config = NULL, verbose = FALSE) {
  # Input validation
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }

  T_obs <- nrow(X)
  N <- ncol(X)

  if (r > N) {
    stop(sprintf("Cannot extract r=%d factors when N=%d predictors available", r, N))
  }
  if (r > T_obs) {
    stop(sprintf("Cannot extract r=%d factors when T=%d observations available", r, T_obs))
  }

  # Optional demeaning (if not already done in preprocessing)
  if (demean) {
    X_centered <- scale(X, center = TRUE, scale = FALSE)
    X_means <- attr(X_centered, "scaled:center")
  } else {
    X_centered <- X
    X_means <- NULL
  }

  if (verbose || isTRUE(config$debug)) {
    log_debug(sprintf(
      "[paper_pca_factors] T=%d, N=%d, r=%d, demean=%s",
      T_obs, N, r, demean
    ), config)
  }

  # Step 1: Compute Σ̂_X = (1/T) X'X
  # This is the sample covariance matrix of predictors (time dimension)
  Sigma_X <- (1 / T_obs) * crossprod(X_centered)  # X'X / T

  # Step 2: Eigen-decompose Σ̂_X
  eigen_result <- eigen(Sigma_X, symmetric = TRUE)
  eigenvalues_all <- eigen_result$values
  eigenvectors_all <- eigen_result$vectors

  # Check for non-positive eigenvalues (should not happen for PSD matrix)
  if (any(eigenvalues_all[1:r] <= 0)) {
    warning("Some eigenvalues are non-positive; covariance matrix may be singular")
  }

  # Select top r eigenvalues and eigenvectors
  eigenvalues_r <- eigenvalues_all[1:r]
  V_r <- eigenvectors_all[, 1:r, drop = FALSE]  # N x r matrix

  # Step 3: Scale eigenvectors to satisfy N^{-1} Λ'Λ = I_r
  # V_r has V_r' V_r = I_r (eigenvectors are orthonormal)
  # We need Λ̂ such that N^{-1} Λ̂' Λ̂ = I_r
  # This means Λ̂' Λ̂ = N I_r
  # So Λ̂ = sqrt(N) * V_r
  Lambda <- sqrt(N) * V_r  # N x r matrix

  # Verify normalization
  norm_check <- (1/N) * crossprod(Lambda)  # Should equal I_r
  I_r <- diag(r)
  max_deviation <- max(abs(norm_check - I_r))

  if (verbose || isTRUE(config$debug)) {
    log_debug(sprintf(
      "[paper_pca_factors] Normalization check: max|N^{-1}Λ'Λ - I| = %.2e (should be ~ 0)",
      max_deviation
    ), config)
  }

  if (max_deviation > 1e-10) {
    warning(sprintf(
      "PCA normalization N^{-1}Λ'Λ deviates from I by %.2e (tolerance 1e-10)",
      max_deviation
    ))
  }

  # Step 4: Compute factors F̂ = N^{-1} X Λ̂
  F_hat <- (1/N) * X_centered %*% Lambda  # T x r matrix

  # Verify reconstruction (optional diagnostic)
  if (verbose || isTRUE(config$debug)) {
    # Reconstruction: X ≈ F̂ Λ̂'
    X_reconstructed <- F_hat %*% t(Lambda)
    reconstruction_error <- sum((X_centered - X_reconstructed)^2) / (N * T_obs)

    # Also check variance explained
    total_variance <- sum(diag(Sigma_X))
    explained_variance <- sum(eigenvalues_r)
    pct_explained <- 100 * explained_variance / total_variance

    log_debug(sprintf(
      "[paper_pca_factors] Reconstruction MSE: %.4f | Variance explained: %.1f%%",
      reconstruction_error, pct_explained
    ), config)
  }

  # Return results
  list(
    F = F_hat,                      # T x r factor scores
    Lambda = Lambda,                # N x r loadings with N^{-1}Λ'Λ = I
    eigenvalues = eigenvalues_r,    # Top r eigenvalues
    Sigma_X = Sigma_X,              # N x N covariance matrix
    method = "paper_pca",           # Identifier
    norm_check_max_dev = max_deviation,  # Diagnostic
    X_means = X_means               # NULL if demean=FALSE
  )
}


#' Wrapper for extract_pca using paper-compliant implementation
#'
#' Drop-in replacement for the current extract_pca() that uses
#' the paper-compliant implementation instead of prcomp().
#'
#' @param X Matrix (T x N) of predictors, already standardized
#' @param k_max Maximum number of principal components to extract
#' @param config Configuration list
#' @param use_paper_method Logical, if TRUE use paper_pca_factors,
#'   if FALSE use original prcomp (default TRUE)
#'
#' @return A list matching extract_pca() interface:
#'   - F: Matrix (T x k_max) of factor scores
#'   - loadings: Matrix (N x k_max) of factor loadings
#'   - pca_obj: The paper PCA object or prcomp object
#'
#' @export
extract_pca_paper <- function(X, k_max = 12, config = NULL, use_paper_method = TRUE) {
  # Input validation
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }

  if (nrow(X) < k_max) {
    stop(sprintf("Not enough observations (%d) to extract %d components", nrow(X), k_max))
  }

  if (ncol(X) < k_max) {
    log_warn(sprintf("Number of variables (%d) < k_max (%d). Using k_max = %d", ncol(X), k_max, ncol(X)), config)
    k_max <- ncol(X)
  }

  log_debug(sprintf("Extracting %d PCA factors from %d x %d matrix (paper method=%s)",
                    k_max, nrow(X), ncol(X), use_paper_method), config)

  if (use_paper_method) {
    # Use paper-compliant implementation
    # Assume X is already centered in preprocessing, so demean=FALSE
    pca_result <- paper_pca_factors(
      X = X,
      r = k_max,
      demean = FALSE,  # X should already be mean-zero from preprocessing
      config = config,
      verbose = isTRUE(config$debug)
    )

    # Return in extract_pca() format
    list(
      F = pca_result$F,
      loadings = pca_result$Lambda,
      pca_obj = pca_result
    )
  } else {
    # Fall back to original prcomp implementation
    pca <- prcomp(X, center = isTRUE(config$pca_center), scale. = isTRUE(config$pca_scale))
    F_hat <- pca$x[, 1:k_max, drop = FALSE]
    loadings <- pca$rotation[, 1:k_max, drop = FALSE]

    log_warn("[extract_pca_paper] Using prcomp (NOT paper-compliant normalization)", config)

    list(
      F = F_hat,
      loadings = loadings,
      pca_obj = pca
    )
  }
}
