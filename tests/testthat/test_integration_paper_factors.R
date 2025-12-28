# Integration Tests: Paper-Compliant Factor Extraction in Forecasting Pipeline
# Tests that the paper implementations work correctly when called through
# extract_factors_at_origin() in the full forecasting workflow

# Source all modules if not already loaded
if (!exists("extract_factors_at_origin")) {
  source_files <- list.files("../../R", pattern = "\\.R$", full.names = TRUE)
  for (f in source_files) source(f, local = FALSE)
}

test_that("Paper PCA integration works in forecasting pipeline", {
  skip_if_not_installed("fbi")
  skip_if_not_installed("pls")

  # Create minimal test configuration
  config <- config_us_default()
  config$use_paper_pca <- TRUE
  config$use_paper_pls <- FALSE  # Test PCA only
  config$debug <- FALSE

  # Create synthetic dataset
  set.seed(123)
  T <- 120
  N <- 50
  dates <- seq.Date(as.Date("2010-01-01"), by = "month", length.out = T)

  # Create panel with predictors
  panel_final <- data.frame(
    date = dates,
    matrix(rnorm(T * N), T, N)
  )
  colnames(panel_final)[-1] <- paste0("X", 1:N)

  # Create target variable
  target_vals <- rnorm(T)
  targets_list <- list(
    TARGET = list(
      h1 = target_vals,
      h6 = c(rep(NA, 5), target_vals[1:(T-5)]),
      h12 = c(rep(NA, 11), target_vals[1:(T-11)])
    )
  )

  # Test extraction at origin t=60
  t_origin <- dates[60]

  result <- extract_factors_at_origin(
    panel_final = panel_final,
    targets_list = targets_list,
    target_name = "TARGET",
    h = 1,
    t_origin = t_origin,
    k_max_pca = 5,
    k_max_pls = 5,
    config = config
  )

  # Verify structure
  expect_true("pca" %in% names(result))
  expect_true("F" %in% names(result$pca))
  expect_true("loadings" %in% names(result$pca))

  # Verify dimensions
  expect_equal(nrow(result$pca$F), 60)  # Up to t_origin
  expect_equal(ncol(result$pca$F), 5)   # k_max_pca = 5
  expect_equal(nrow(result$pca$loadings), N)
  expect_equal(ncol(result$pca$loadings), 5)

  # Verify paper normalization: N^{-1}Λ'Λ = I_r
  Lambda <- result$pca$loadings
  norm_check <- (1/N) * crossprod(Lambda)
  I_r <- diag(5)
  max_dev <- max(abs(norm_check - I_r))
  expect_lt(max_dev, 1e-10, label = "PCA normalization N^{-1}Λ'Λ = I_r")
})

test_that("Paper PLS integration works in forecasting pipeline", {
  skip_if_not_installed("fbi")
  skip_if_not_installed("pls")

  # Create minimal test configuration
  config <- config_us_default()
  config$use_paper_pca <- FALSE  # Test PLS only
  config$use_paper_pls <- TRUE
  config$debug <- FALSE

  # Create synthetic dataset
  set.seed(456)
  T <- 120
  N <- 50
  dates <- seq.Date(as.Date("2010-01-01"), by = "month", length.out = T)

  panel_final <- data.frame(
    date = dates,
    matrix(rnorm(T * N), T, N)
  )
  colnames(panel_final)[-1] <- paste0("X", 1:N)

  # Create correlated target variable
  X_mat <- as.matrix(panel_final[, -1])
  true_factor <- X_mat %*% rnorm(N)
  target_vals <- true_factor + rnorm(T, sd = 0.5)

  targets_list <- list(
    TARGET = list(
      h1 = target_vals,
      h6 = c(rep(NA, 5), target_vals[1:(T-5)])
    )
  )

  # Test extraction at origin t=60
  t_origin <- dates[60]

  result <- extract_factors_at_origin(
    panel_final = panel_final,
    targets_list = targets_list,
    target_name = "TARGET",
    h = 1,
    t_origin = t_origin,
    k_max_pca = 5,
    k_max_pls = 5,
    config = config
  )

  # Verify structure
  expect_true("pls" %in% names(result))
  expect_true("F" %in% names(result$pls))

  # Verify dimensions
  expect_equal(nrow(result$pls$F), 60)
  expect_equal(ncol(result$pls$F), 5)

  # Verify paper normalization: N^{-1}α'α = 1
  # Extract weights from the model
  pls_model <- result$pls$model
  expect_true("weights" %in% names(pls_model))

  weights <- pls_model$weights
  for (j in 1:5) {
    alpha_j <- weights[, j]
    constraint_val <- (1/N) * sum(alpha_j^2)
    expect_equal(constraint_val, 1.0, tolerance = 1e-10,
                 label = sprintf("PLS constraint N^{-1}α'α = 1 for factor %d", j))
  }

  # Verify orthogonality of factors
  F_mat <- result$pls$F
  F_cross <- crossprod(F_mat) / 60
  off_diag <- F_cross - diag(diag(F_cross))
  max_off_diag <- max(abs(off_diag))
  expect_lt(max_off_diag, 1e-6, label = "PLS factors are approximately orthogonal")
})

test_that("Both paper methods work together in pipeline", {
  skip_if_not_installed("fbi")
  skip_if_not_installed("pls")

  # Create minimal test configuration
  config <- config_us_default()
  config$use_paper_pca <- TRUE
  config$use_paper_pls <- TRUE
  config$debug <- TRUE  # Enable to test logging

  # Create synthetic dataset
  set.seed(789)
  T <- 100
  N <- 40
  dates <- seq.Date(as.Date("2010-01-01"), by = "month", length.out = T)

  panel_final <- data.frame(
    date = dates,
    matrix(rnorm(T * N), T, N)
  )
  colnames(panel_final)[-1] <- paste0("X", 1:N)

  target_vals <- rnorm(T)
  targets_list <- list(
    TARGET = list(h1 = target_vals)
  )

  # Test extraction
  t_origin <- dates[60]

  result <- extract_factors_at_origin(
    panel_final = panel_final,
    targets_list = targets_list,
    target_name = "TARGET",
    h = 1,
    t_origin = t_origin,
    k_max_pca = 8,
    k_max_pls = 8,
    config = config
  )

  # Both PCA and PLS should be present
  expect_true(all(c("pca", "pls") %in% names(result)))

  # Both should have correct dimensions
  expect_equal(ncol(result$pca$F), 8)
  expect_equal(ncol(result$pls$F), 8)
  expect_equal(nrow(result$pca$F), nrow(result$pls$F))
})

test_that("Fallback to library methods works when paper flags are FALSE", {
  skip_if_not_installed("fbi")
  skip_if_not_installed("pls")

  # Create configuration with library methods
  config <- config_us_default()
  config$use_paper_pca <- FALSE
  config$use_paper_pls <- FALSE
  config$debug <- FALSE

  # Create synthetic dataset
  set.seed(999)
  T <- 80
  N <- 30
  dates <- seq.Date(as.Date("2010-01-01"), by = "month", length.out = T)

  panel_final <- data.frame(
    date = dates,
    matrix(rnorm(T * N), T, N)
  )
  colnames(panel_final)[-1] <- paste0("X", 1:N)

  target_vals <- rnorm(T)
  targets_list <- list(
    TARGET = list(h1 = target_vals)
  )

  t_origin <- dates[50]

  # Should not error when using library methods
  expect_no_error({
    result <- extract_factors_at_origin(
      panel_final = panel_final,
      targets_list = targets_list,
      target_name = "TARGET",
      h = 1,
      t_origin = t_origin,
      k_max_pca = 5,
      k_max_pls = 5,
      config = config
    )
  })

  # Structure should still be correct
  expect_true(all(c("pca", "pls") %in% names(result)))
  expect_equal(ncol(result$pca$F), 5)
  expect_equal(ncol(result$pls$F), 5)
})

test_that("Paper implementations handle horizon-specific targets correctly", {
  skip_if_not_installed("fbi")
  skip_if_not_installed("pls")

  config <- config_us_default()
  config$use_paper_pca <- TRUE
  config$use_paper_pls <- TRUE
  config$debug <- FALSE

  set.seed(111)
  T <- 120
  N <- 40
  dates <- seq.Date(as.Date("2010-01-01"), by = "month", length.out = T)

  panel_final <- data.frame(
    date = dates,
    matrix(rnorm(T * N), T, N)
  )
  colnames(panel_final)[-1] <- paste0("X", 1:N)

  target_base <- rnorm(T)

  # Create horizon-specific targets (as done in actual preprocessing)
  targets_list <- list(
    TARGET = list(
      h1 = target_base,
      h6 = c(rep(NA, 5), target_base[1:(T-5)]),
      h12 = c(rep(NA, 11), target_base[1:(T-11)]),
      h24 = c(rep(NA, 23), target_base[1:(T-23)])
    )
  )

  t_origin <- dates[60]

  # Test different horizons
  for (h in c(1, 6, 12, 24)) {
    result <- extract_factors_at_origin(
      panel_final = panel_final,
      targets_list = targets_list,
      target_name = "TARGET",
      h = h,
      t_origin = t_origin,
      k_max_pca = 4,
      k_max_pls = 4,
      config = config
    )

    # PLS should extract factors using horizon-specific y
    expect_equal(ncol(result$pls$F), 4,
                 label = sprintf("PLS should extract 4 factors at h=%d", h))

    # PCA should be independent of horizon
    expect_equal(ncol(result$pca$F), 4,
                 label = sprintf("PCA should extract 4 factors at h=%d", h))
  }
})
