# Comprehensive Tests for Paper-Compliant Factor Extraction
# Tests verify exact compliance with Bae (2024) Section 2.2

library(testthat)

# Source the paper implementations
# (These should be loaded via package mechanism, but explicitly source for testing)
source("../../R/paper_factors_pca.R")
source("../../R/paper_factors_pls.R")
source("../../R/utils_logging.R")

# =============================================================================
# PCA Tests: Verify Bae (2024) Section 2.2 Compliance
# =============================================================================

test_that("Paper PCA: Normalization N^{-1}Λ'Λ = I_r is satisfied", {
  set.seed(123)
  T_obs <- 100
  N <- 50
  r <- 5

  X <- matrix(rnorm(T_obs * N), T_obs, N)
  X <- scale(X, center = TRUE, scale = FALSE)  # Mean-zero

  pca_result <- paper_pca_factors(X, r = r, demean = FALSE)

  Lambda <- pca_result$Lambda
  I_r <- diag(r)
  norm_matrix <- (1/N) * crossprod(Lambda)  # N^{-1} Λ'Λ

  # Should be identity within numerical tolerance
  expect_equal(norm_matrix, I_r, tolerance = 1e-10)

  # Explicit check
  max_dev <- max(abs(norm_matrix - I_r))
  expect_lt(max_dev, 1e-10)
})

test_that("Paper PCA: Factor formula F̂ = N^{-1} X Λ̂ is correct", {
  set.seed(456)
  T_obs <- 80
  N <- 40
  r <- 3

  X <- matrix(rnorm(T_obs * N), T_obs, N)
  X <- scale(X, center = TRUE, scale = FALSE)

  pca_result <- paper_pca_factors(X, r = r, demean = FALSE)

  # Verify F̂ = N^{-1} X Λ̂
  F_computed <- (1/N) * X %*% pca_result$Lambda
  F_returned <- pca_result$F

  expect_equal(F_computed, F_returned, tolerance = 1e-12)
})

test_that("Paper PCA: Reconstruction X ≈ F̂ Λ̂' holds", {
  set.seed(789)
  T_obs <- 100
  N <- 30
  r <- 10

  X <- matrix(rnorm(T_obs * N), T_obs, N)
  X <- scale(X, center = TRUE, scale = FALSE)

  pca_result <- paper_pca_factors(X, r = r, demean = FALSE)

  # Reconstruction
  X_reconstructed <- pca_result$F %*% t(pca_result$Lambda)

  # Should not be exactly equal (dimension reduction), but should be good approximation
  # Check Frobenius norm of residual
  residual <- X - X_reconstructed
  frobenius_norm <- sqrt(sum(residual^2))

  # Reconstruction should explain substantial variance
  total_variance <- sum(X^2)
  reconstruction_variance <- sum(X_reconstructed^2)

  expect_gt(reconstruction_variance / total_variance, 0.5)  # At least 50% variance explained
})

test_that("Paper PCA: Eigenvalues are in descending order", {
  set.seed(111)
  T_obs <- 60
  N <- 25
  r <- 5

  X <- matrix(rnorm(T_obs * N), T_obs, N)
  X <- scale(X, center = TRUE, scale = FALSE)

  pca_result <- paper_pca_factors(X, r = r, demean = FALSE)

  eigenvalues <- pca_result$eigenvalues

  # Should be in descending order
  expect_true(all(diff(eigenvalues) <= 0))

  # Should all be positive (PSD matrix)
  expect_true(all(eigenvalues > 0))
})

test_that("Paper PCA: Covariance matrix Σ̂_X = (1/T) X'X is correct", {
  set.seed(222)
  T_obs <- 70
  N <- 20
  r <- 3

  X <- matrix(rnorm(T_obs * N), T_obs, N)
  X <- scale(X, center = TRUE, scale = FALSE)

  pca_result <- paper_pca_factors(X, r = r, demean = FALSE)

  Sigma_X_computed <- (1/T_obs) * crossprod(X)
  Sigma_X_returned <- pca_result$Sigma_X

  expect_equal(Sigma_X_computed, Sigma_X_returned, tolerance = 1e-14)
})

test_that("Paper PCA vs prcomp: Different normalizations confirmed", {
  set.seed(333)
  T_obs <- 50
  N <- 30
  r <- 4

  X <- matrix(rnorm(T_obs * N), T_obs, N)
  X <- scale(X, center = TRUE, scale = FALSE)

  # Paper PCA
  pca_paper <- paper_pca_factors(X, r = r, demean = FALSE)
  Lambda_paper <- pca_paper$Lambda

  # Base R prcomp
  pca_base <- prcomp(X, center = FALSE, scale. = FALSE)
  V_base <- pca_base$rotation[, 1:r]

  # V_base should satisfy V'V = I
  V_cross <- crossprod(V_base)
  expect_equal(V_cross, diag(r), tolerance = 1e-12)

  # Lambda_paper should satisfy N^{-1}Λ'Λ = I
  Lambda_cross <- (1/N) * crossprod(Lambda_paper)
  expect_equal(Lambda_cross, diag(r), tolerance = 1e-12)

  # Relationship: Lambda_paper ≈ sqrt(N) * V_base (up to sign flips)
  for (j in 1:r) {
    # Allow for sign flip
    ratio1 <- Lambda_paper[1, j] / V_base[1, j]
    ratio2 <- -Lambda_paper[1, j] / V_base[1, j]

    # One of these should be approximately sqrt(N)
    expect_true(abs(ratio1 - sqrt(N)) < 0.01 || abs(ratio2 - sqrt(N)) < 0.01)
  }
})

# =============================================================================
# PLS Tests: Verify Bae (2024) Table 1 Compliance
# =============================================================================

test_that("Paper PLS: Normalization N^{-1}α'α = 1 for all factors", {
  set.seed(444)
  T_obs <- 100
  N <- 40
  k <- 5

  X <- matrix(rnorm(T_obs * N), T_obs, N)
  y <- rnorm(T_obs)
  X <- scale(X, center = TRUE, scale = FALSE)
  y <- y - mean(y)

  pls_result <- paper_pls_factors(X, y, k = k, demean_X = FALSE, demean_y = FALSE)

  weights <- pls_result$weights

  for (j in 1:k) {
    alpha_j <- weights[, j]
    constraint_val <- (1/N) * sum(alpha_j^2)

    # Should equal 1 within tolerance
    expect_equal(constraint_val, 1.0, tolerance = 1e-10)
  }
})

test_that("Paper PLS: Deflation creates orthogonal factors", {
  set.seed(555)
  T_obs <- 120
  N <- 30
  k <- 6

  X <- matrix(rnorm(T_obs * N), T_obs, N)
  y <- rnorm(T_obs)
  X <- scale(X, center = TRUE, scale = FALSE)
  y <- y - mean(y)

  pls_result <- paper_pls_factors(X, y, k = k, demean_X = FALSE, demean_y = FALSE)

  F_matrix <- pls_result$F

  # Compute F'F / T
  F_cross <- crossprod(F_matrix) / T_obs

  # Diagonal elements (variances) should be positive
  expect_true(all(diag(F_cross) > 0))

  # Off-diagonal elements (correlations * std devs) should be near zero
  off_diag <- F_cross - diag(diag(F_cross))
  max_off_diag <- max(abs(off_diag))

  # Allow small numerical error
  expect_lt(max_off_diag, 1e-10)
})

test_that("Paper PLS: First factor matches F̂_1 = X α1 formula", {
  set.seed(666)
  T_obs <- 80
  N <- 25
  k <- 1  # Just first factor

  X <- matrix(rnorm(T_obs * N), T_obs, N)
  y <- rnorm(T_obs)
  X <- scale(X, center = TRUE, scale = FALSE)
  y <- y - mean(y)

  pls_result <- paper_pls_factors(X, y, k = k, demean_X = FALSE, demean_y = FALSE)

  # First factor
  F_1 <- pls_result$F[, 1]
  alpha_1 <- pls_result$weights[, 1]

  # Should satisfy F̂_1 = X α1
  F_1_computed <- X %*% alpha_1

  expect_equal(F_1_computed, F_1, tolerance = 1e-12)

  # Also check α1 ∝ X'y with correct normalization
  v <- crossprod(X, y)  # X'y
  v_norm <- sqrt(sum(v^2))
  alpha_1_expected <- sqrt(N) * (v / v_norm)

  expect_equal(alpha_1, as.vector(alpha_1_expected), tolerance = 1e-12)
})

test_that("Paper PLS: Deflated X is orthogonal to previous factors", {
  set.seed(777)
  T_obs <- 100
  N <- 30
  k <- 4

  X <- matrix(rnorm(T_obs * N), T_obs, N)
  y <- rnorm(T_obs)
  X <- scale(X, center = TRUE, scale = FALSE)
  y <- y - mean(y)

  pls_result <- paper_pls_factors(X, y, k = k, demean_X = FALSE, demean_y = FALSE, verbose = FALSE)

  # Check orthogonality residuals from diagnostics
  ortho_checks <- pls_result$diagnostics$orthogonality_residuals

  # Each should be near zero (F_prev' X* ≈ 0)
  for (j in 1:(k-1)) {
    expect_lt(ortho_checks[[j]], 1e-10)
  }
})

test_that("Paper PLS: F̂_j is in column space of X*_j (deflated X)", {
  set.seed(888)
  T_obs <- 90
  N <- 20
  k <- 3

  X <- matrix(rnorm(T_obs * N), T_obs, N)
  y <- rnorm(T_obs)
  X <- scale(X, center = TRUE, scale = FALSE)
  y <- y - mean(y)

  pls_result <- paper_pls_factors(X, y, k = k, demean_X = FALSE, demean_y = FALSE)

  # For each factor j, verify F̂_j = X*_j α_j
  for (j in 1:k) {
    X_star_j <- pls_result$X_deflated_seq[[j]]
    alpha_j <- pls_result$weights[, j]
    F_j_expected <- X_star_j %*% alpha_j
    F_j_actual <- pls_result$F[, j]

    expect_equal(F_j_expected, F_j_actual, tolerance = 1e-12)
  }
})

test_that("Paper PLS: Diagnostics report correct deviations", {
  set.seed(999)
  T_obs <- 110
  N <- 35
  k <- 5

  X <- matrix(rnorm(T_obs * N), T_obs, N)
  y <- rnorm(T_obs)
  X <- scale(X, center = TRUE, scale = FALSE)
  y <- y - mean(y)

  pls_result <- paper_pls_factors(X, y, k = k, demean_X = FALSE, demean_y = FALSE)

  diagnostics <- pls_result$diagnostics

  # Constraint deviations should all be near zero
  expect_true(all(diagnostics$constraint_deviations < 1e-10))

  # Final off-diagonal max should be near zero
  expect_lt(diagnostics$final_off_diagonal_max, 1e-10)

  # F cross product should be approximately diagonal
  F_cross <- diagnostics$F_cross_product
  diag_vals <- diag(F_cross)
  expect_true(all(diag_vals > 0))  # Positive variances

  off_diag <- F_cross - diag(diag_vals)
  expect_lt(max(abs(off_diag)), 1e-10)
})

# =============================================================================
# Integration Tests: Wrappers
# =============================================================================

test_that("extract_pca_paper wrapper works correctly", {
  skip_if_not(exists("extract_pca_paper"), "extract_pca_paper not loaded")

  set.seed(1010)
  T_obs <- 80
  N <- 40
  k_max <- 8

  X <- matrix(rnorm(T_obs * N), T_obs, N)
  X <- scale(X, center = TRUE, scale = FALSE)

  config <- list(pca_center = FALSE, pca_scale = FALSE, debug = FALSE)

  result <- extract_pca_paper(X, k_max = k_max, config = config, use_paper_method = TRUE)

  # Check structure
  expect_true(is.list(result))
  expect_true("F" %in% names(result))
  expect_true("loadings" %in% names(result))
  expect_true("pca_obj" %in% names(result))

  # Check dimensions
  expect_equal(nrow(result$F), T_obs)
  expect_equal(ncol(result$F), k_max)
  expect_equal(nrow(result$loadings), N)
  expect_equal(ncol(result$loadings), k_max)

  # Verify paper normalization
  Lambda <- result$loadings
  norm_matrix <- (1/N) * crossprod(Lambda)
  expect_equal(norm_matrix, diag(k_max), tolerance = 1e-10)
})

test_that("extract_pls_paper wrapper works correctly", {
  skip_if_not(exists("extract_pls_paper"), "extract_pls_paper not loaded")

  set.seed(1111)
  T_obs <- 100
  N <- 30
  k_max <- 6

  X <- matrix(rnorm(T_obs * N), T_obs, N)
  y <- rnorm(T_obs)
  X <- scale(X, center = TRUE, scale = FALSE)
  y <- y - mean(y)

  config <- list(pls_center = FALSE, pls_scale = FALSE, pls_method = "oscorespls", debug = FALSE)

  result <- extract_pls_paper(X, y, k_max = k_max, config = config, use_paper_method = TRUE)

  # Check structure
  expect_true(is.list(result))
  expect_true("F" %in% names(result))
  expect_true("model" %in% names(result))

  # Check dimensions
  expect_equal(nrow(result$F), T_obs)
  expect_equal(ncol(result$F), k_max)

  # Verify constraint on weights
  weights <- result$model$weights
  for (j in 1:k_max) {
    alpha_j <- weights[, j]
    constraint_val <- (1/N) * sum(alpha_j^2)
    expect_equal(constraint_val, 1.0, tolerance = 1e-10)
  }
})

test_that("Paper PLS handles NA values in y correctly", {
  skip_if_not(exists("extract_pls_paper"), "extract_pls_paper not loaded")

  set.seed(1212)
  T_obs <- 100
  N <- 25
  k_max <- 4

  X <- matrix(rnorm(T_obs * N), T_obs, N)
  y <- rnorm(T_obs)
  X <- scale(X, center = TRUE, scale = FALSE)
  y <- y - mean(y, na.rm = TRUE)

  # Introduce some NAs in y
  y[c(95:100)] <- NA

  config <- list(pls_center = FALSE, pls_scale = FALSE, debug = FALSE)

  # Should not error; should fit on non-NA observations
  expect_error({
    result <- extract_pls_paper(X, y, k_max = k_max, config = config, use_paper_method = TRUE)
  }, NA)

  result <- extract_pls_paper(X, y, k_max = k_max, config = config, use_paper_method = TRUE)

  # Should still return factors for all T observations
  expect_equal(nrow(result$F), T_obs)
})

# =============================================================================
# Comparison Tests: Paper vs Library Implementations
# =============================================================================

test_that("PCA: Paper implementation differs from prcomp as expected", {
  set.seed(1313)
  T_obs <- 60
  N <- 20
  r <- 5

  X <- matrix(rnorm(T_obs * N), T_obs, N)
  X <- scale(X, center = TRUE, scale = FALSE)

  # Paper PCA
  pca_paper <- paper_pca_factors(X, r = r, demean = FALSE)

  # prcomp PCA
  pca_base <- prcomp(X, center = FALSE, scale. = FALSE)
  F_base <- pca_base$x[, 1:r]
  Lambda_base <- pca_base$rotation[, 1:r]

  # Normalizations should differ
  # prcomp: Λ'Λ = I
  # paper: N^{-1}Λ'Λ = I

  Lambda_base_norm <- crossprod(Lambda_base)
  Lambda_paper_norm <- (1/N) * crossprod(pca_paper$Lambda)

  expect_equal(Lambda_base_norm, diag(r), tolerance = 1e-12)
  expect_equal(Lambda_paper_norm, diag(r), tolerance = 1e-12)

  # But absolute scales differ
  expect_false(isTRUE(all.equal(Lambda_base, pca_paper$Lambda, tolerance = 1e-6)))

  # Factors also differ in scale
  expect_false(isTRUE(all.equal(F_base, pca_paper$F, tolerance = 1e-6)))
})

test_that("PLS: Paper implementation can be compared to pls::plsr scores", {
  skip_if_not_installed("pls")
  set.seed(1414)
  T_obs <- 100
  N <- 30
  k <- 3

  X <- matrix(rnorm(T_obs * N), T_obs, N)
  y <- rnorm(T_obs)
  X <- scale(X, center = TRUE, scale = FALSE)
  y <- y - mean(y)

  # Paper PLS
  pls_paper <- paper_pls_factors(X, y, k = k, demean_X = FALSE, demean_y = FALSE)

  # Library PLS
  df <- data.frame(y = y, X)
  pls_lib <- pls::plsr(y ~ ., data = df, ncomp = k, center = FALSE, scale = FALSE, method = "oscorespls")
  scores_lib <- pls_lib$scores

  # They may differ in scale/sign, but should span similar spaces
  # This is a weaker test - just check both produce reasonable results

  expect_equal(nrow(pls_paper$F), T_obs)
  expect_equal(ncol(pls_paper$F), k)
  expect_equal(nrow(scores_lib), T_obs)
  expect_equal(ncol(scores_lib), k)

  # Both should have non-trivial variance
  expect_true(all(apply(pls_paper$F, 2, sd) > 0.01))
  expect_true(all(apply(scores_lib, 2, sd) > 0.01))
})

# =============================================================================
# Edge Cases and Error Handling
# =============================================================================

test_that("Paper PCA handles r > N error", {
  set.seed(1515)
  T_obs <- 50
  N <- 10
  r <- 15  # More than N

  X <- matrix(rnorm(T_obs * N), T_obs, N)

  expect_error(
    paper_pca_factors(X, r = r),
    "Cannot extract r=15 factors when N=10"
  )
})

test_that("Paper PLS handles k > N error", {
  set.seed(1616)
  T_obs <- 60
  N <- 10
  k <- 15

  X <- matrix(rnorm(T_obs * N), T_obs, N)
  y <- rnorm(T_obs)

  expect_error(
    paper_pls_factors(X, y, k = k),
    "Cannot extract k=15 PLS factors when N=10"
  )
})

test_that("Paper PLS handles low-variance y gracefully", {
  set.seed(1717)
  T_obs <- 80
  N <- 20
  k <- 5

  X <- matrix(rnorm(T_obs * N), T_obs, N)
  y <- rep(0, T_obs)  # Zero variance

  # Should still run but may produce warnings
  expect_warning({
    result <- paper_pls_factors(X, y, k = k, demean_X = FALSE, demean_y = FALSE)
  }, "covariance exhausted")
})

# =============================================================================
# Summary message
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("Paper Factor Extraction Tests Summary:\n")
cat("================================================================================\n")
cat("All tests verify exact compliance with Bae (2024) Section 2.2 specifications:\n")
cat("  - PCA: N^{-1}Λ'Λ = I_r normalization\n")
cat("  - PCA: F̂ = N^{-1} X Λ̂ factor formula\n")
cat("  - PCA: Σ̂_X = (1/T) X'X covariance\n")
cat("  - PLS: N^{-1}α'α = 1 constraint\n")
cat("  - PLS: Q(F_prev)X deflation procedure\n")
cat("  - PLS: Orthogonality of factors\n")
cat("================================================================================\n")
cat("\n")
