# Tests for k-selection rules (BN-BIC and Onatski)

test_that("compute_bn_bic_k returns valid structure", {
  set.seed(42)
  T_obs <- 100
  N <- 50
  X <- matrix(rnorm(T_obs * N), T_obs, N)

  result <- compute_bn_bic_k(X, k_max = 10, min_k = 1L)

  expect_type(result, "list")
  expect_true("k_hat" %in% names(result))
  expect_true("BIC_values" %in% names(result))
  expect_true("V_values" %in% names(result))
  expect_true(result$k_hat >= 1)
  expect_true(result$k_hat <= 10)
  expect_length(result$BIC_values, 10)
  expect_length(result$V_values, 10)
})

test_that("compute_bn_bic_k enforces min_k = 1", {
  set.seed(123)
  # Pure noise - BN-BIC might want k=0, but we enforce min_k=1
  T_obs <- 100
  N <- 20
  X <- matrix(rnorm(T_obs * N), T_obs, N)

  result <- compute_bn_bic_k(X, k_max = 5, min_k = 1L)

  expect_true(result$k_hat >= 1)
  # If k_hat_raw was 0, we should have a warning
  if (result$k_hat_raw < 1) {
    expect_false(is.null(result$warning))
  }
})

test_that("compute_bn_bic_k identifies factors in synthetic data", {
  set.seed(42)
  T_obs <- 200
  N <- 50
  k_true <- 3

  # Generate factor model: X = F * Lambda' + E
  F_true <- matrix(rnorm(T_obs * k_true), T_obs, k_true)
  Lambda <- matrix(rnorm(N * k_true), N, k_true)
  E <- matrix(rnorm(T_obs * N) * 0.5, T_obs, N)  # Small noise
  X <- F_true %*% t(Lambda) + E

  result <- compute_bn_bic_k(X, k_max = 10, min_k = 1L)

  # Should identify approximately k_true factors (allowing some tolerance)
  expect_true(result$k_hat >= k_true - 1)
  expect_true(result$k_hat <= k_true + 2)
})

test_that("compute_bn_bic_k handles edge cases", {
  set.seed(42)

  # k_max = 1
  X <- matrix(rnorm(100 * 20), 100, 20)
  result <- compute_bn_bic_k(X, k_max = 1, min_k = 1L)
  expect_equal(result$k_hat, 1L)

  # Small T (T < N)
  X_small_T <- matrix(rnorm(30 * 50), 30, 50)
  result_small <- compute_bn_bic_k(X_small_T, k_max = 10, min_k = 1L)
  expect_true(result_small$k_hat >= 1)
})

test_that("compute_onatski_k returns valid structure", {
  set.seed(42)
  T_obs <- 100
  N <- 50
  X <- matrix(rnorm(T_obs * N), T_obs, N)

  result <- compute_onatski_k(X, r_max = 10, min_k = 1L)

  expect_type(result, "list")
  expect_true("k_hat" %in% names(result))
  expect_true("eigenvalues" %in% names(result))
  expect_true("u_hat" %in% names(result))
  expect_true("delta" %in% names(result))
  expect_true(result$k_hat >= 1)
})

test_that("compute_onatski_k enforces min_k = 1", {
  set.seed(123)
  # Pure noise - Onatski might want k=0
  T_obs <- 100
  N <- 50
  X <- matrix(rnorm(T_obs * N), T_obs, N)

  result <- compute_onatski_k(X, r_max = 10, min_k = 1L)

  expect_true(result$k_hat >= 1)
})

test_that("compute_onatski_k auto-adjusts r_max", {
  set.seed(42)
  # Small dimensions where r_max = 12 would violate constraint
  T_obs <- 20
  N <- 15
  X <- matrix(rnorm(T_obs * N), T_obs, N)

  # r_max = 12 should be auto-adjusted since 2*12+1 > min(15, 20)
  result <- compute_onatski_k(X, r_max = 12, min_k = 1L)

  expect_true(result$r_max_used < 12)
  expect_true(result$r_max_used >= 1)
  expect_true(result$k_hat >= 1)
})

test_that("compute_onatski_k identifies factors in synthetic data", {
  set.seed(42)
  T_obs <- 200
  N <- 50
  k_true <- 2

  # Generate factor model with clear eigenvalue gap
  F_true <- matrix(rnorm(T_obs * k_true) * 3, T_obs, k_true)  # Strong factors
  Lambda <- matrix(rnorm(N * k_true), N, k_true)
  E <- matrix(rnorm(T_obs * N) * 0.3, T_obs, N)  # Small noise
  X <- F_true %*% t(Lambda) + E

  result <- compute_onatski_k(X, r_max = 10, min_k = 1L)

  # Should identify approximately k_true factors
  expect_true(result$k_hat >= 1)
  expect_true(result$k_hat <= k_true + 2)
})

test_that("compute_cov_eigenvalues returns correct dimensions", {
  set.seed(42)
  T_obs <- 100
  N <- 30
  X <- matrix(rnorm(T_obs * N), T_obs, N)

  eigenvalues <- compute_cov_eigenvalues(X)

  expect_length(eigenvalues, N)
  expect_true(all(eigenvalues >= 0))  # Eigenvalues should be non-negative
  expect_true(all(diff(eigenvalues) <= 0))
})

test_that("select_k_dynamic dispatches correctly", {
  set.seed(42)
  T_obs <- 100
  N <- 50
  X <- matrix(rnorm(T_obs * N), T_obs, N)

  config <- list(
    k_selection_settings = list(
      min_k = 1,
      bn_bic_sigma_sq = "v_kmax",
      onatski_r_max = 10
    )
  )

  # BN-BIC
  result_bn <- select_k_dynamic(X, "bn_bic", k_max = 8, config = config)
  expect_true(result_bn$k_hat >= 1)
  expect_equal(result_bn$k_rule, "bn_bic")

  # Onatski
  result_on <- select_k_dynamic(X, "onatski", k_max = 8, config = config)
  expect_true(result_on$k_hat >= 1)
  expect_equal(result_on$k_rule, "onatski")
})

test_that("k_hat is always <= k_max", {
  set.seed(42)
  T_obs <- 100
  N <- 50

  # Strong factor model that might suggest many factors
  k_true <- 8
  F_true <- matrix(rnorm(T_obs * k_true) * 5, T_obs, k_true)
  Lambda <- matrix(rnorm(N * k_true), N, k_true)
  E <- matrix(rnorm(T_obs * N) * 0.1, T_obs, N)
  X <- F_true %*% t(Lambda) + E

  # Test with k_max = 4 (less than true k)
  result_bn <- compute_bn_bic_k(X, k_max = 4, min_k = 1L)
  expect_true(result_bn$k_hat <= 4)

  # Onatski with k_max = 4 should clamp result
  result_on <- compute_onatski_k(X, r_max = 10, min_k = 1L, k_max = 4L)
  expect_true(result_on$k_hat >= 1)
  expect_true(result_on$k_hat <= 4)

  # Onatski without k_max (for comparison)
  result_on_unclamped <- compute_onatski_k(X, r_max = 10, min_k = 1L, k_max = NULL)
  expect_true(result_on_unclamped$k_hat >= 1)
})

test_that("BN-BIC V_values are decreasing", {
  set.seed(42)
  T_obs <- 100
  N <- 50
  X <- matrix(rnorm(T_obs * N), T_obs, N)

  result <- compute_bn_bic_k(X, k_max = 10, min_k = 1L)

  # V(k) should be decreasing in k (more factors = less residual variance)
  expect_true(all(diff(result$V_values) <= 0))
})
