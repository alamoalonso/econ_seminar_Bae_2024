# Smoke tests for basic functionality
# These tests verify that the core workflow runs without errors

library(testthat)

test_that("Configuration system works", {
  config <- config_us_default()

  expect_true(is.list(config))
  expect_equal(config$dataset_id, "US_FRED")
  expect_true(all(c("horizons", "k_max_pca", "k_max_pls") %in% names(config)))

  # Validation should pass
  expect_silent(validate_config(config))
})

test_that("Configuration update works", {
  config <- config_us_default()
  config_updated <- update_config(config, debug = TRUE, horizons = c(1))

  expect_true(config_updated$debug)
  expect_equal(config_updated$horizons, c(1))
})

test_that("Outlier removal works", {
  x <- c(1, 2, 3, 4, 5, 100)  # 100 is an outlier
  x_clean <- remove_outliers_iqr(x, multiplier = 10)

  expect_true(is.na(x_clean[6]))
  expect_equal(x_clean[1:5], c(1, 2, 3, 4, 5))
})

test_that("Target construction works", {
  x <- 1:10
  h <- 2

  y_h <- construct_target_h(x, h)

  expect_equal(length(y_h), 10)
  expect_equal(y_h[1], 3)   # x[1+2]
  expect_equal(y_h[2], 4)   # x[2+2]
  expect_true(is.na(y_h[9]))  # Beyond end
  expect_true(is.na(y_h[10])) # Beyond end
})

test_that("PCA extraction works", {
  set.seed(123)
  X <- matrix(rnorm(100*20), 100, 20)

  pca_result <- extract_pca(X, k_max = 5, config = list(pca_center = TRUE, pca_scale = FALSE))

  expect_equal(nrow(pca_result$F), 100)
  expect_equal(ncol(pca_result$F), 5)
  expect_equal(nrow(pca_result$loadings), 20)
  expect_equal(ncol(pca_result$loadings), 5)
})

test_that("PLS extraction works", {
  set.seed(123)
  X <- matrix(rnorm(100*20), 100, 20)
  y <- rnorm(100)

  pls_result <- extract_pls(X, y, k_max = 5,
                            config = list(pls_center = TRUE, pls_scale = FALSE, pls_method = "oscorespls"))

  expect_equal(nrow(pls_result$F), 100)
  expect_equal(ncol(pls_result$F), 5)
})

test_that("AR model fitting and prediction works", {
  set.seed(123)
  n <- 100
  y_base <- arima.sim(list(ar = c(0.5, 0.3)), n = n)
  y_h <- c(y_base[3:n], NA, NA)  # h=2 target

  sample_idx <- 1:80
  fit <- fit_AR_BIC(y_h, y_base, sample_idx, max_p = 6)

  expect_false(is.null(fit))
  expect_true("model" %in% names(fit))
  expect_true("p" %in% names(fit))

  pred <- predict_AR_at_origin(fit, y_base, t_index = 81)
  expect_true(is.finite(pred))
})

test_that("Window functions work", {
  dates <- seq.Date(as.Date("2000-01-01"), by = "month", length.out = 240)

  # Recursive
  idx_rec <- get_recursive_idx(120)
  expect_equal(idx_rec, 1:120)

  # Rolling (10 years = 120 months)
  idx_roll <- get_rolling_idx(dates, 240)
  expect_true(length(idx_roll) <= 240)
  expect_true(all(idx_roll <= 240))
})

test_that("BIC computation works", {
  residuals <- c(0.1, -0.2, 0.15, -0.1, 0.05)
  bic <- compute_bic(residuals, n = 5, k_param = 3)

  expect_true(is.finite(bic))
  expect_true(bic > 0)  # BIC should be positive for typical cases
})

# Integration test (requires current.csv to exist)
test_that("Full workflow runs (if data available)", {
  skip_if_not(file.exists("current.csv"), "current.csv not found, skipping integration test")

  # Create minimal config
  config <- config_us_default()
  config$series_list <- c("INDPRO")  # Single series
  config$horizons <- c(1)             # Single horizon
  config$factor_methods <- c("PCA")   # Single method
  config$debug <- FALSE

  # This should run without errors
  expect_error({
    results <- run_workflow(config)
    rmse_results <- compute_rmse(results, config)
  }, NA)  # NA means "expect no error"
})
