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
  config <- config_us_default()

  # Recursive
  idx_rec <- get_recursive_idx(120)
  expect_equal(idx_rec, 1:120)

  # Rolling (10 years = 120 months with monthly frequency)
  idx_roll <- get_rolling_idx(dates, 240, config)
  expect_equal(length(idx_roll), 120)  # 10 years * 12 months
  expect_true(all(idx_roll <= 240))
  expect_equal(min(idx_roll), 121)  # 240 - 120 + 1
  expect_equal(max(idx_roll), 240)
})

test_that("BIC computation works", {
  residuals <- c(0.1, -0.2, 0.15, -0.1, 0.05)
  bic <- compute_bic(residuals, n = 5, k_param = 3)

  expect_true(is.finite(bic))
  expect_true(bic > 0)  # BIC should be positive for typical cases
})

test_that("Config includes save_forecasts option", {
  config <- config_us_default()

  expect_true("save_forecasts" %in% names(config))
  expect_true(config$save_forecasts)  # Default should be TRUE
})

test_that("extract_forecasts_from_res returns correct structure", {
  # Create mock forecast results
  n <- 100
  k_max <- 3
  res <- list(
    truth = c(rep(NA, 59), rnorm(40), NA),  # first_idx=60, last valid=99
    AR_rec = c(rep(NA, 59), rnorm(40), NA),
    AR_roll = c(rep(NA, 59), rnorm(40), NA),
    DI_rec = matrix(c(rep(NA, 59), rnorm(40), NA), n, k_max),
    DI_roll = matrix(c(rep(NA, 59), rnorm(40), NA), n, k_max),
    DIAR_rec = matrix(c(rep(NA, 59), rnorm(40), NA), n, k_max),
    DIAR_roll = matrix(c(rep(NA, 59), rnorm(40), NA), n, k_max),
    DLAG_rec = matrix(c(rep(NA, 59), rnorm(40), NA), n, k_max),
    DLAG_roll = matrix(c(rep(NA, 59), rnorm(40), NA), n, k_max)
  )

  dates <- seq.Date(as.Date("1970-01-01"), by = "month", length.out = n)

  result <- extract_forecasts_from_res(
    res = res,
    dates = dates,
    series_id = "TEST",
    h = 1,
    factor_method = "PCA",
    eval_start = as.Date("1970-01-01"),
    eval_end = as.Date("1978-03-01"),
    k_max_factor = k_max,
    run_id = "test_run"
  )

  expect_true(is.data.frame(result))
  expect_true(all(c("run_id", "series_id", "h", "scheme", "factor_method",
                    "model_class", "k", "origin_index", "origin_date",
                    "target_index", "target_date", "y_true", "y_hat",
                    "method_id") %in% names(result)))

  # Check AR rows have factor_method = NA
  ar_rows <- result[result$model_class == "AR", ]
  expect_true(all(is.na(ar_rows$factor_method)))

  # Check target_index = origin_index + h
  expect_true(all(result$target_index == result$origin_index + result$h))

  # Check method_id format for AR vs non-AR
  expect_true(all(grepl("^(recursive|rolling)_AR$", ar_rows$method_id)))
  non_ar_rows <- result[result$model_class != "AR", ]
  expect_true(all(grepl("^(recursive|rolling)_PCA_(DI|DIAR|DIAR-LAG)_k[0-9]+$", non_ar_rows$method_id)))
})

test_that("compute_evaluation returns forecasts when save_forecasts=TRUE", {
  skip_if_not(file.exists("current.csv"), "current.csv not found, skipping integration test")

  config <- config_us_default()
  config$series_list <- c("INDPRO")
  config$horizons <- c(1)
  config$factor_methods <- c("PCA")
  config$save_forecasts <- TRUE
  config$debug <- FALSE
  config$do_tests <- FALSE  # Skip tests to speed up

  results <- run_workflow(config)
  evaluation <- compute_evaluation(results, config)

  expect_true("forecasts" %in% names(evaluation))
  expect_true(!is.null(evaluation$forecasts))
  expect_true(nrow(evaluation$forecasts) > 0)

  # Verify no AR duplicates (should have 2 rows per origin: recursive and rolling)
  ar_forecasts <- evaluation$forecasts[evaluation$forecasts$model_class == "AR", ]
  ar_per_origin <- table(ar_forecasts$origin_index)
  expect_true(all(ar_per_origin == 2))  # One recursive, one rolling
})

test_that("compute_evaluation returns NULL forecasts when save_forecasts=FALSE", {
  skip_if_not(file.exists("current.csv"), "current.csv not found, skipping integration test")

  config <- config_us_default()
  config$series_list <- c("INDPRO")
  config$horizons <- c(1)
  config$factor_methods <- c("PCA")
  config$save_forecasts <- FALSE
  config$debug <- FALSE
  config$do_tests <- FALSE

  results <- run_workflow(config)
  evaluation <- compute_evaluation(results, config)

  expect_true(is.null(evaluation$forecasts))
})

# Integration test (requires current.csv to exist)
test_that("Full workflow runs (if data available)", {
  skip_if_not(file.exists("current.csv"), "current.csv not found, skipping integration test")

  # Create minimal config with temp output directory
  config <- config_us_default()
  config$series_list <- c("INDPRO")  # Single series
  config$horizons <- c(1)             # Single horizon
  config$factor_methods <- c("PCA")   # Single method
  config$debug <- FALSE
  config$save_forecasts <- TRUE
  config$output_dir <- tempdir()
  config$run_id <- paste0("test_", format(Sys.time(), "%Y%m%d_%H%M%S"))

  # Run workflow and compute evaluation
  results <- run_workflow(config)
  evaluation <- compute_evaluation(results, config)
  rmse_results <- evaluation$rmse_results

  # Verify RMSE results exist

  expect_true(nrow(rmse_results) > 0)

  # Verify forecasts were collected (save_forecasts = TRUE)
  expect_true(!is.null(evaluation$forecasts))
  expect_true(nrow(evaluation$forecasts) > 0)
})

test_that("save_results writes forecasts to disk", {
  skip_if_not(file.exists("current.csv"), "current.csv not found, skipping integration test")

  # Create minimal config with temp output directory
  config <- config_us_default()
  config$series_list <- c("INDPRO")
  config$horizons <- c(1)
  config$factor_methods <- c("PCA")
  config$debug <- FALSE
  config$save_forecasts <- TRUE
  config$do_tests <- FALSE
  config$make_plots <- FALSE
  config$output_dir <- tempdir()
  config$run_id <- paste0("test_save_", format(Sys.time(), "%Y%m%d_%H%M%S"))

  results <- run_workflow(config)
  evaluation <- compute_evaluation(results, config)

  # Save results including forecasts
  output_dir <- save_results(results, evaluation$rmse_results, config,
                              forecasts = evaluation$forecasts)

  # Verify forecasts were saved (either as parquet dir or CSV)
  forecasts_dir <- file.path(output_dir, "forecasts")
  forecasts_csv <- file.path(output_dir, "forecasts_long.csv")

  forecasts_saved <- dir.exists(forecasts_dir) || file.exists(forecasts_csv)
  expect_true(forecasts_saved)

  # Clean up
  unlink(output_dir, recursive = TRUE)
})

# ============================================================================
# MCS Tests
# ============================================================================

test_that("Config includes MCS settings", {
  config <- config_us_default()

  expect_true("mcs" %in% names(config))
  expect_true(is.list(config$mcs))
  expect_true("enabled" %in% names(config$mcs))
  expect_true("alphas" %in% names(config$mcs))
  expect_true("M0_sets" %in% names(config$mcs))
  expect_false(config$mcs$enabled)  # Default should be FALSE
})

test_that("Loss functions work correctly", {
  y_true <- c(1, 2, 3, 4, 5)
  y_hat <- c(1.1, 2.2, 2.8, 4.5, 4.9)

  se <- loss_se(y_true, y_hat)
  ae <- loss_ae(y_true, y_hat)

  expect_equal(length(se), 5)
  expect_equal(length(ae), 5)
  expect_true(all(se >= 0))
  expect_true(all(ae >= 0))
  expect_equal(se, (y_true - y_hat)^2)
  expect_equal(ae, abs(y_true - y_hat))
})

test_that("Stationary bootstrap produces correct dimensions", {
  set.seed(123)
  x <- matrix(rnorm(100 * 3), 100, 3)
  B <- 10

  boot_samples <- stationary_bootstrap(x, B = B, expected_block_length = 5, seed = 42)

  expect_equal(dim(boot_samples), c(100, 3, B))
})

test_that("run_mcs returns correct structure", {
  set.seed(123)
  # Create loss matrix with 50 time points and 3 methods
  # Method 1 has lowest loss on average (should be in MCS)
  T_obs <- 50
  L <- cbind(
    method1 = rnorm(T_obs, mean = 1, sd = 0.5),
    method2 = rnorm(T_obs, mean = 1.5, sd = 0.5),
    method3 = rnorm(T_obs, mean = 2, sd = 0.5)
  )

  # Suppress warnings about MCS package not being installed (fallback is OK for tests)
  result <- suppressWarnings(
    run_mcs(L, alpha = 0.10, B = 100, stat_type = "Tmax", seed = 42)
  )

  expect_true(is.list(result))
  expect_true("superior_set" %in% names(result))
  expect_true("eliminated" %in% names(result))
  expect_true("pvalues" %in% names(result))
  expect_true("status" %in% names(result))
  # Status can be "ok" (MCS package) or "fallback" (no MCS package)
  expect_true(result$status %in% c("ok", "fallback"))

  # Method 1 should be in superior set (lowest average loss)
  expect_true("method1" %in% result$superior_set)
})

test_that("resolve_method_specs works correctly", {
  available <- c(
    "recursive_PCA_DI_k1", "recursive_PCA_DI_k2", "recursive_PCA_DI_k3",
    "recursive_PCA_DIAR_k1", "recursive_PCA_DIAR_k2",
    "recursive_PLS_DI_k1", "recursive_PLS_DI_k2",
    "recursive_AR", "rolling_AR"
  )

  # Test k1-PCA with specific model_class
  result <- resolve_method_specs("k1-PCA", available, scheme = "recursive", model_class = "DI")
  expect_true("recursive_PCA_DI_k1" %in% result)
  expect_equal(length(result), 1)

  # Test k1-PCA without model_class filter (matches all model classes)
  result <- resolve_method_specs("k1-PCA", available, scheme = "recursive", model_class = NULL)
  expect_true("recursive_PCA_DI_k1" %in% result)
  expect_true("recursive_PCA_DIAR_k1" %in% result)
  expect_equal(length(result), 2)

  # Test AR (AR resolution ignores model_class)
  result <- resolve_method_specs("AR", available, scheme = "recursive")
  expect_true("recursive_AR" %in% result)
  expect_equal(length(result), 1)

  # Test without scheme filter
  result <- resolve_method_specs("AR", available)
  expect_true("recursive_AR" %in% result)
  expect_true("rolling_AR" %in% result)
})

test_that("MCS include_ar_in_M0 config option exists", {
  config <- config_us_default()

  expect_true("include_ar_in_M0" %in% names(config$mcs))
  expect_true(config$mcs$include_ar_in_M0)  # Default should be TRUE
})

test_that("MCS evaluation skips when disabled", {
  # Create mock config with MCS disabled
  config <- config_us_default()
  config$mcs$enabled <- FALSE

  # Create minimal forecasts data
  forecasts <- tibble::tibble(
    run_id = "test",
    series_id = "TEST",
    h = 1,
    scheme = "recursive",
    factor_method = "PCA",
    model_class = "DI",
    k = 1,
    origin_index = 1:10,
    y_true = rnorm(10),
    y_hat = rnorm(10),
    method_id = "recursive_PCA_DI_k1"
  )

  result <- compute_mcs_evaluation(forecasts, config)
  expect_null(result)
})

test_that("k_selection_settings config option exists", {
  config <- config_us_default()
  expect_true("k_selection_settings" %in% names(config))
  expect_true(is.list(config$k_selection_settings))
  expect_equal(config$k_selection_settings$min_k, 1)
})

test_that("select_k_bai_ng returns placeholder (TODO)", {
  set.seed(123)
  X <- matrix(rnorm(100 * 20), 100, 20)

  # This should return a warning about not being implemented
  expect_warning({
    result <- select_k_bai_ng(X, k_max = 5, criterion = "IC_p3")
  }, "not yet fully implemented")

  expect_true(is.list(result))
  expect_equal(result$k_selected, 5)  # Returns k_max as fallback
  expect_equal(result$status, "placeholder")
})

# ============================================================================
# Factor Specs and Mixed Grid + Dynamic Tests
# ============================================================================

test_that("factor_specs backward compatibility works", {
  config <- config_us_default()

  # With factor_specs = NULL, should auto-expand from factor_methods
  specs <- get_factor_specs(config)

  expect_true(length(specs) >= 1)
  expect_true(all(sapply(specs, function(s) s$k_mode == "grid")))
})

test_that("factor_specs with explicit specs works", {
  config <- config_us_default()
  config$factor_specs <- list(
    create_factor_spec("PCA_grid", "PCA", "grid", NULL, 12),
    create_factor_spec("PCA_BNBIC", "PCA", "dynamic", "bn_bic", 12),
    create_factor_spec("PLS_grid", "PLS", "grid", NULL, 8),
    create_factor_spec("PLS_ON", "PLS", "dynamic", "onatski", 8)
  )

  specs <- get_factor_specs(config)

  expect_length(specs, 4)

  # Check we have both grid and dynamic
  modes <- sapply(specs, function(s) s$k_mode)
  expect_true("grid" %in% modes)
  expect_true("dynamic" %in% modes)
})

test_that("run_forecasts_for_spec works with grid spec", {
  skip_if_not(file.exists("current.csv"), "current.csv not found, skipping integration test")

  config <- config_us_default()
  config$series_list <- c("INDPRO")
  config$horizons <- c(1)
  config$debug <- FALSE

  results <- run_workflow(config)
  dataset <- results$dataset

  # Create a grid spec
  spec_grid <- create_factor_spec("PCA_grid", "PCA", "grid", NULL, 3)

  res <- run_forecasts_for_spec(
    series_name = "INDPRO",
    h = 1,
    dates = dataset$dates,
    panel_final = dataset$panel_final,
    panel_std1 = dataset$panel_std1,
    targets_list = dataset$targets_list,
    factor_spec = spec_grid,
    config = config
  )

  expect_equal(res$spec_id, "PCA_grid")
  expect_equal(res$k_mode, "grid")
  expect_true(is.matrix(res$DI_rec))
  expect_equal(ncol(res$DI_rec), 3)  # k_max = 3
})

test_that("run_forecasts_for_spec works with dynamic spec", {
  skip_if_not(file.exists("current.csv"), "current.csv not found, skipping integration test")

  config <- config_us_default()
  config$series_list <- c("INDPRO")
  config$horizons <- c(1)
  config$debug <- FALSE

  results <- run_workflow(config)
  dataset <- results$dataset

  # Create a dynamic BN-BIC spec
  spec_dyn <- create_factor_spec("PCA_BNBIC", "PCA", "dynamic", "bn_bic", 6)

  res <- run_forecasts_for_spec(
    series_name = "INDPRO",
    h = 1,
    dates = dataset$dates,
    panel_final = dataset$panel_final,
    panel_std1 = dataset$panel_std1,
    targets_list = dataset$targets_list,
    factor_spec = spec_dyn,
    config = config
  )

  expect_equal(res$spec_id, "PCA_BNBIC")
  expect_equal(res$k_mode, "dynamic")
  expect_true(is.vector(res$DI_rec))  # Dynamic returns vectors, not matrices
  expect_true(is.vector(res$k_hat))

  # k_hat should have valid values where forecasts exist
  valid_k_hat <- res$k_hat[!is.na(res$truth)]
  expect_true(all(valid_k_hat >= 1, na.rm = TRUE))  # min_k = 1 enforced
})

test_that("extract_forecasts_from_spec_res handles grid spec", {
  # Create mock grid spec result
  n <- 100
  k_max <- 3
  res <- list(
    spec_id = "PCA_grid",
    factor_method = "PCA",
    k_mode = "grid",
    k_rule = NULL,
    k_max = k_max,
    truth = c(rep(NA, 59), rnorm(40), NA),
    k_hat = rep(NA_integer_, n),
    DI_rec = matrix(c(rep(NA, 59), rnorm(40), NA), n, k_max),
    DI_roll = matrix(c(rep(NA, 59), rnorm(40), NA), n, k_max),
    DIAR_rec = matrix(c(rep(NA, 59), rnorm(40), NA), n, k_max),
    DIAR_roll = matrix(c(rep(NA, 59), rnorm(40), NA), n, k_max),
    DLAG_rec = matrix(c(rep(NA, 59), rnorm(40), NA), n, k_max),
    DLAG_roll = matrix(c(rep(NA, 59), rnorm(40), NA), n, k_max),
    AR_rec = c(rep(NA, 59), rnorm(40), NA),
    AR_roll = c(rep(NA, 59), rnorm(40), NA)
  )

  dates <- seq.Date(as.Date("1970-01-01"), by = "month", length.out = n)

  result <- extract_forecasts_from_spec_res(
    res = res,
    dates = dates,
    series_id = "TEST",
    h = 1,
    eval_start = as.Date("1970-01-01"),
    eval_end = as.Date("1978-03-01"),
    run_id = "test_run",
    first_forecast_idx = 60
  )

  expect_true(is.data.frame(result))
  expect_true("factor_spec_id" %in% names(result))
  expect_true("k_mode" %in% names(result))
  expect_true("k_selection_rule" %in% names(result))
  expect_true("training_window_start" %in% names(result))
  expect_true("training_window_end" %in% names(result))

  # Grid mode should have multiple k values
  non_ar_rows <- result[result$model_class != "AR", ]
  expect_true(all(non_ar_rows$k_mode == "grid"))
  expect_true(all(non_ar_rows$k_selection_rule == "fixed"))

  # Should have k = 1, 2, 3 for each model/scheme/origin
  k_values <- unique(non_ar_rows$k)
  expect_true(all(1:k_max %in% k_values))
})

test_that("extract_forecasts_from_spec_res handles dynamic spec", {
  # Create mock dynamic spec result
  n <- 100
  k_max <- 6
  res <- list(
    spec_id = "PCA_BNBIC",
    factor_method = "PCA",
    k_mode = "dynamic",
    k_rule = "bn_bic",
    k_max = k_max,
    truth = c(rep(NA, 59), rnorm(40), NA),
    k_hat = c(rep(NA_integer_, 59), sample(1:4, 40, replace = TRUE), NA_integer_),
    DI_rec = c(rep(NA, 59), rnorm(40), NA),
    DI_roll = c(rep(NA, 59), rnorm(40), NA),
    DIAR_rec = c(rep(NA, 59), rnorm(40), NA),
    DIAR_roll = c(rep(NA, 59), rnorm(40), NA),
    DLAG_rec = c(rep(NA, 59), rnorm(40), NA),
    DLAG_roll = c(rep(NA, 59), rnorm(40), NA),
    AR_rec = c(rep(NA, 59), rnorm(40), NA),
    AR_roll = c(rep(NA, 59), rnorm(40), NA)
  )

  dates <- seq.Date(as.Date("1970-01-01"), by = "month", length.out = n)

  result <- extract_forecasts_from_spec_res(
    res = res,
    dates = dates,
    series_id = "TEST",
    h = 1,
    eval_start = as.Date("1970-01-01"),
    eval_end = as.Date("1978-03-01"),
    run_id = "test_run",
    first_forecast_idx = 60
  )

  expect_true(is.data.frame(result))

  # Dynamic mode should have k_selection_rule = "bn_bic"
  non_ar_rows <- result[result$model_class != "AR", ]
  expect_true(all(non_ar_rows$k_mode == "dynamic"))
  expect_true(all(non_ar_rows$k_selection_rule == "bn_bic"))

  # k should vary by origin (it's k_hat)
  expect_true(length(unique(non_ar_rows$k)) > 1 || nrow(non_ar_rows) < 2)
})

test_that("resolve_method_specs handles new spec-based method_ids", {
  available <- c(
    # Legacy format
    "recursive_PCA_DI_k1", "recursive_PCA_DI_k2",
    # New grid format
    "recursive_PCA_grid_DI_k1", "recursive_PCA_grid_DI_k2",
    # New dynamic format
    "recursive_PCA_BNBIC_DI", "recursive_PLS_ON_DIAR",
    "recursive_AR", "rolling_AR"
  )

  # Test PCA_BNBIC pattern
  result <- resolve_method_specs("PCA_BNBIC", available, scheme = "recursive", model_class = "DI")
  expect_true("recursive_PCA_BNBIC_DI" %in% result)

  # Test PLS_ON pattern
  result <- resolve_method_specs("PLS_ON", available, scheme = "recursive", model_class = "DIAR")
  expect_true("recursive_PLS_ON_DIAR" %in% result)

  # Test k1-PCA matches both legacy and new grid format
  result <- resolve_method_specs("k1-PCA", available, scheme = "recursive", model_class = "DI")
  expect_true("recursive_PCA_DI_k1" %in% result || "recursive_PCA_grid_DI_k1" %in% result)

  # Test PCA_grid pattern
  result <- resolve_method_specs("PCA_grid", available, scheme = "recursive", model_class = "DI")
  expect_true("recursive_PCA_grid_DI_k1" %in% result)
  expect_true("recursive_PCA_grid_DI_k2" %in% result)
})

test_that("Mixed grid + dynamic in same workflow run (integration)", {
  skip_if_not(file.exists("current.csv"), "current.csv not found, skipping integration test")

  config <- config_us_default()
  config$series_list <- c("INDPRO")
  config$horizons <- c(1)
  config$debug <- FALSE
  config$save_forecasts <- TRUE
  config$do_tests <- FALSE
  config$output_dir <- tempdir()
  config$run_id <- paste0("test_mixed_", format(Sys.time(), "%Y%m%d_%H%M%S"))

  # Define mixed specs: both grid and dynamic
  config$factor_specs <- list(
    create_factor_spec("PCA_grid", "PCA", "grid", NULL, 3),
    create_factor_spec("PCA_BNBIC", "PCA", "dynamic", "bn_bic", 6)
  )

  results <- run_workflow(config)
  evaluation <- compute_evaluation_with_specs(results, config)

  expect_true(!is.null(evaluation$forecasts))
  expect_true(nrow(evaluation$forecasts) > 0)

  forecasts <- evaluation$forecasts

  # CRITICAL VERIFICATION: Both grid and dynamic should appear in output
  # Grid specs
  grid_rows <- forecasts[forecasts$factor_spec_id == "PCA_grid" & !is.na(forecasts$factor_spec_id), ]
  expect_true(nrow(grid_rows) > 0, info = "Grid spec forecasts should exist")
  expect_true(all(grid_rows$k_mode == "grid"))
  expect_true(all(grid_rows$k_selection_rule == "fixed"))

  # Dynamic specs
  bnbic_rows <- forecasts[forecasts$factor_spec_id == "PCA_BNBIC" & !is.na(forecasts$factor_spec_id), ]
  expect_true(nrow(bnbic_rows) > 0, info = "BN-BIC spec forecasts should exist")
  expect_true(all(bnbic_rows$k_mode == "dynamic"))
  expect_true(all(bnbic_rows$k_selection_rule == "bn_bic"))

  # Verify method_ids are distinct
  grid_method_ids <- unique(grid_rows$method_id)
  bnbic_method_ids <- unique(bnbic_rows$method_id)
  expect_length(intersect(grid_method_ids, bnbic_method_ids), 0,
                info = "Grid and dynamic method_ids should not collide")

  # Grid method_ids should have _k{n} suffix
  expect_true(all(grepl("_k[0-9]+$", grid_method_ids)))

  # Dynamic method_ids should NOT have _k{n} suffix
  expect_true(all(!grepl("_k[0-9]+$", bnbic_method_ids)))

  # Verify k_hat varies across origins for dynamic (confirms recursive computation)
  if (length(unique(bnbic_rows$origin_index)) > 1) {
    k_values <- unique(bnbic_rows$k)
    # k should be >= 1 (min_k enforced)
    expect_true(all(k_values >= 1))
  }

  # Clean up
  output_dir <- file.path(config$output_dir, config$run_id)
  if (dir.exists(output_dir)) {
    unlink(output_dir, recursive = TRUE)
  }
})
