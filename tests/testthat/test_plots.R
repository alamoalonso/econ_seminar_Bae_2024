# Tests for plotting functions
# These tests verify that the RMSE plotting functions work correctly

library(testthat)

test_that("summarise_mean_rmse_by_k works with synthetic data", {
  # Create synthetic RMSE table
  rmse_tbl <- data.frame(
    series = rep(c("A", "B", "C"), each = 48),
    h = rep(rep(c(1, 6, 12, 24), each = 12), 3),
    scheme = rep(rep(c("recursive", "rolling"), each = 6), 12),
    factor_method = rep(c("PCA", "PLS"), 72),
    model = rep(c("DI", "DIAR", "DIAR-LAG"), each = 4, times = 12),
    k = rep(c(1, 2, 3, 4), times = 36),
    mse = runif(144, 0.5, 2),
    mse_ar = runif(144, 0.8, 1.5),
    rmse_rel = runif(144, 0.8, 1.2)
  )

  # Test summarization
  summary <- summarise_mean_rmse_by_k(
    rmse_tbl,
    metric = "rmse_rel",
    models = c("DI", "DIAR", "DIAR-LAG"),
    horizons = c(1, 6, 12, 24)
  )

  expect_true(is.data.frame(summary))
  expect_true(all(c("scheme", "factor_method", "model", "h", "k", "mean_rmse", "n_series") %in% names(summary)))
  expect_true(all(summary$n_series > 0))
  expect_true(all(is.finite(summary$mean_rmse)))
})

test_that("DIAR-LAG k constraint is enforced", {
  # Create synthetic RMSE table with k > 4 for DIAR-LAG
  rmse_tbl <- data.frame(
    series = rep("A", 10),
    h = rep(1, 10),
    scheme = rep("recursive", 10),
    factor_method = rep("PCA", 10),
    model = rep("DIAR-LAG", 10),
    k = 1:10,  # k goes up to 10
    mse = runif(10, 0.5, 2),
    mse_ar = 1,
    rmse_rel = runif(10, 0.8, 1.2)
  )

  # Test that only k <= 4 is retained
  summary <- summarise_mean_rmse_by_k(
    rmse_tbl,
    metric = "rmse_rel",
    models = c("DIAR-LAG"),
    horizons = c(1),
    k_max_diarlag = 4
  )

  expect_true(max(summary$k) <= 4)
  expect_equal(nrow(summary), 4)  # Should only have k=1,2,3,4
})

test_that("plot_mean_rmse_kfactor returns patchwork object", {
  # Create synthetic RMSE table with varying scales (DI much larger than others)
  rmse_tbl <- data.frame(
    series = rep(c("A", "B"), each = 96),
    h = rep(rep(c(1, 6, 12, 24), each = 24), 2),
    scheme = rep(rep(c("recursive", "rolling"), each = 12), 8),
    factor_method = rep("PCA", 192),
    model = rep(rep(c("DI", "DIAR", "DIAR-LAG"), each = 4), 16),
    k = rep(rep(c(1, 2, 3, 4), times = 3), 16),
    mse = runif(192, 0.5, 2),
    mse_ar = 1,
    rmse_rel = runif(192, 0.8, 1.2)
  )

  # Make DI have larger values to test independent scaling
  rmse_tbl$rmse_rel[rmse_tbl$model == "DI"] <- runif(sum(rmse_tbl$model == "DI"), 2, 5)

  # Test plot creation
  p <- plot_mean_rmse_kfactor(
    rmse_tbl,
    factor_method = "PCA",
    metric = "rmse_rel",
    schemes = c("recursive", "rolling"),
    models = c("DI", "DIAR", "DIAR-LAG"),
    horizons = c(1, 6, 12, 24),
    save_path = NULL
  )

  # Should return a patchwork object (combination of ggplots)
  expect_true(inherits(p, c("patchwork", "gg", "ggplot")))
})

test_that("plot_mean_rmse_kfactor handles missing data gracefully", {
  # Create synthetic RMSE table with only PLS data
  rmse_tbl <- data.frame(
    series = rep("A", 10),
    h = rep(1, 10),
    scheme = rep("recursive", 10),
    factor_method = rep("PLS", 10),  # Only PLS, no PCA
    model = rep("DI", 10),
    k = 1:10,
    mse = runif(10, 0.5, 2),
    mse_ar = 1,
    rmse_rel = runif(10, 0.8, 1.2)
  )

  # Test that requesting PCA gives warning and returns NULL
  expect_warning({
    p <- plot_mean_rmse_kfactor(
      rmse_tbl,
      factor_method = "PCA",
      metric = "rmse_rel",
      save_path = NULL
    )
  }, "No data available")

  expect_null(suppressWarnings(
    plot_mean_rmse_kfactor(rmse_tbl, factor_method = "PCA", save_path = NULL)
  ))
})

test_that("plot_bae_fig1_kpca creates file", {
  skip_on_ci <- function() {
    skip_if(Sys.getenv("CI") != "", "Skipping on CI")
  }
  skip_on_ci()

  # Create synthetic RMSE table
  rmse_tbl <- data.frame(
    series = rep(c("A", "B"), each = 96),
    h = rep(rep(c(1, 6, 12, 24), each = 24), 2),
    scheme = rep(rep(c("recursive", "rolling"), each = 12), 8),
    factor_method = rep("PCA", 192),
    model = rep(rep(c("DI", "DIAR", "DIAR-LAG"), each = 4), 16),
    k = rep(rep(c(1, 2, 3, 4), times = 3), 16),
    mse = runif(192, 0.5, 2),
    mse_ar = 1,
    rmse_rel = runif(192, 0.8, 1.2)
  )

  # Create temporary directory
  temp_dir <- tempdir()

  # Test file creation
  plot_bae_fig1_kpca(
    rmse_tbl,
    metric = "rmse_rel",
    save_dir = temp_dir
  )

  expect_true(file.exists(file.path(temp_dir, "fig1_mean_rmse_k_pca.png")))

  # Clean up
  unlink(file.path(temp_dir, "fig1_mean_rmse_k_pca.png"))
})

test_that("plot_bae_fig2_kpls creates file", {
  skip_on_ci <- function() {
    skip_if(Sys.getenv("CI") != "", "Skipping on CI")
  }
  skip_on_ci()

  # Create synthetic RMSE table
  rmse_tbl <- data.frame(
    series = rep(c("A", "B"), each = 96),
    h = rep(rep(c(1, 6, 12, 24), each = 24), 2),
    scheme = rep(rep(c("recursive", "rolling"), each = 12), 8),
    factor_method = rep("PLS", 192),
    model = rep(rep(c("DI", "DIAR", "DIAR-LAG"), each = 4), 16),
    k = rep(rep(c(1, 2, 3, 4), times = 3), 16),
    mse = runif(192, 0.5, 2),
    mse_ar = 1,
    rmse_rel = runif(192, 0.8, 1.2)
  )

  # Create temporary directory
  temp_dir <- tempdir()

  # Test file creation
  plot_bae_fig2_kpls(
    rmse_tbl,
    metric = "rmse_rel",
    save_dir = temp_dir
  )

  expect_true(file.exists(file.path(temp_dir, "fig2_mean_rmse_k_pls.png")))

  # Clean up
  unlink(file.path(temp_dir, "fig2_mean_rmse_k_pls.png"))
})

test_that("summarise_mean_rmse_by_k validates required columns", {
  # Create incomplete table missing rmse_rel
  rmse_tbl <- data.frame(
    series = "A",
    h = 1,
    scheme = "recursive",
    factor_method = "PCA",
    model = "DI",
    k = 1
  )

  expect_error(
    summarise_mean_rmse_by_k(rmse_tbl),
    "Missing required columns"
  )
})

test_that("generate_plots_for_run works with saved results", {
  skip_on_ci <- function() {
    skip_if(Sys.getenv("CI") != "", "Skipping on CI")
  }
  skip_on_ci()

  # Create temporary directory structure
  temp_dir <- tempdir()
  run_dir <- file.path(temp_dir, "test_run_20231214")
  dir.create(run_dir, showWarnings = FALSE)

  # Create synthetic RMSE results
  rmse_tbl <- data.frame(
    series = rep(c("A", "B"), each = 96),
    h = rep(rep(c(1, 6, 12, 24), each = 24), 2),
    scheme = rep(rep(c("recursive", "rolling"), each = 12), 8),
    factor_method = rep(c("PCA", "PLS"), each = 96),
    model = rep(rep(c("DI", "DIAR", "DIAR-LAG"), each = 4), 16),
    k = rep(rep(c(1, 2, 3, 4), times = 3), 16),
    mse = runif(192, 0.5, 2),
    mse_ar = 1,
    rmse_rel = runif(192, 0.8, 1.2)
  )

  # Save to CSV
  write.csv(rmse_tbl, file.path(run_dir, "rmse_results.csv"), row.names = FALSE)

  # Test plot generation
  result <- generate_plots_for_run(run_dir)

  expect_true(file.exists(file.path(run_dir, "fig1_mean_rmse_k_pca.png")))
  expect_true(file.exists(file.path(run_dir, "fig2_mean_rmse_k_pls.png")))
  expect_true(is.list(result))
  expect_equal(length(result), 2)  # Both PCA and PLS plots

  # Clean up
  unlink(run_dir, recursive = TRUE)
})

test_that("generate_plots_for_run handles missing directory gracefully", {
  expect_error(
    generate_plots_for_run("nonexistent_directory"),
    "Run directory not found"
  )
})

test_that("generate_plots_for_run handles missing rmse_results.csv", {
  temp_dir <- tempdir()
  run_dir <- file.path(temp_dir, "empty_run")
  dir.create(run_dir, showWarnings = FALSE)

  expect_error(
    generate_plots_for_run(run_dir),
    "rmse_results.csv not found"
  )

  # Clean up
  unlink(run_dir, recursive = TRUE)
})
