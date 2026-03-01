# Tests for factor_specs module

test_that("create_factor_spec validates inputs correctly", {
  # Valid grid spec
  spec <- create_factor_spec("PCA_grid", "PCA", "grid", NULL, 12)
  expect_equal(spec$id, "PCA_grid")
  expect_equal(spec$factor_method, "PCA")
  expect_equal(spec$k_mode, "grid")
  expect_null(spec$k_rule)
  expect_equal(spec$k_max, 12L)

  # Valid dynamic spec
  spec_dyn <- create_factor_spec("PCA_BNBIC", "PCA", "dynamic", "bn_bic", 12)
  expect_equal(spec_dyn$k_rule, "bn_bic")
  expect_equal(spec_dyn$k_mode, "dynamic")

  # Invalid factor_method

  expect_error(
    create_factor_spec("test", "INVALID", "grid", NULL, 12),
    "Invalid factor_method"
  )

  # Invalid k_mode
  expect_error(
    create_factor_spec("test", "PCA", "invalid", NULL, 12),
    "Invalid k_mode"
  )

  # k_rule provided for grid mode
  expect_error(
    create_factor_spec("test", "PCA", "grid", "bn_bic", 12),
    "k_rule must be NULL for k_mode = 'grid'"
  )

  # Missing k_rule for dynamic mode
  expect_error(
    create_factor_spec("test", "PCA", "dynamic", NULL, 12),
    "k_rule must be one of"
  )

  # Wrong rule for factor_method
  expect_error(
    create_factor_spec("test", "PLS", "dynamic", "bn_bic", 12),
    "bn_bic rule is only valid for PCA"
  )
  expect_error(
    create_factor_spec("test", "PCA", "dynamic", "onatski", 12),
    "onatski rule is only valid for PLS"
  )

  # Invalid k_max
  expect_error(
    create_factor_spec("test", "PCA", "grid", NULL, 0),
    "k_max must be a positive integer"
  )
  expect_error(
    create_factor_spec("test", "PCA", "grid", NULL, -1),
    "k_max must be a positive integer"
  )
})

test_that("expand_factor_methods_to_specs creates correct grid specs", {
  specs <- expand_factor_methods_to_specs(c("PCA", "PLS"), k_max_pca = 12, k_max_pls = 8)

  expect_length(specs, 2)
  expect_equal(names(specs), c("PCA_grid", "PLS_grid"))

  expect_equal(specs$PCA_grid$factor_method, "PCA")
  expect_equal(specs$PCA_grid$k_mode, "grid")
  expect_equal(specs$PCA_grid$k_max, 12)

  expect_equal(specs$PLS_grid$factor_method, "PLS")
  expect_equal(specs$PLS_grid$k_mode, "grid")
  expect_equal(specs$PLS_grid$k_max, 8)
})

test_that("validate_factor_specs catches errors", {
  # Empty list
  expect_error(
    validate_factor_specs(list()),
    "non-empty list"
  )

  # Missing fields
  expect_error(
    validate_factor_specs(list(list(id = "test"))),
    "missing required fields"
  )

  # Duplicate ids
  specs <- list(
    list(id = "dup", factor_method = "PCA", k_mode = "grid", k_rule = NULL, k_max = 12),
    list(id = "dup", factor_method = "PLS", k_mode = "grid", k_rule = NULL, k_max = 12)
  )
  expect_error(
    validate_factor_specs(specs),
    "Duplicate factor_spec ids"
  )
})

test_that("get_factor_specs handles backward compatibility", {
  # Config with factor_methods only (backward compat)
  config_old <- list(
    factor_methods = c("PCA", "PLS"),
    k_max_pca = 10,
    k_max_pls = 8,
    factor_specs = NULL
  )

  specs <- get_factor_specs(config_old)
  expect_length(specs, 2)
  expect_equal(specs[[1]]$k_mode, "grid")

  # Config with explicit factor_specs
  config_new <- list(
    factor_specs = list(
      list(id = "PCA_grid", factor_method = "PCA", k_mode = "grid", k_rule = NULL, k_max = 12),
      list(id = "PCA_BNBIC", factor_method = "PCA", k_mode = "dynamic", k_rule = "bn_bic", k_max = 12)
    )
  )

  specs_new <- get_factor_specs(config_new)
  expect_length(specs_new, 2)
  expect_equal(specs_new[[1]]$k_mode, "grid")
  expect_equal(specs_new[[2]]$k_mode, "dynamic")
})

test_that("get_max_k_for_method works correctly", {
  specs <- list(
    list(id = "PCA_grid", factor_method = "PCA", k_mode = "grid", k_rule = NULL, k_max = 12),
    list(id = "PCA_BNBIC", factor_method = "PCA", k_mode = "dynamic", k_rule = "bn_bic", k_max = 8),
    list(id = "PLS_grid", factor_method = "PLS", k_mode = "grid", k_rule = NULL, k_max = 10)
  )

  expect_equal(get_max_k_for_method(specs, "PCA"), 12)  # max(12, 8)
  expect_equal(get_max_k_for_method(specs, "PLS"), 10)
  expect_equal(get_max_k_for_method(specs, "1-PLS"), 0)  # No matching specs
})

test_that("has_dynamic_rule works correctly", {
  specs_mixed <- list(
    list(id = "PCA_grid", factor_method = "PCA", k_mode = "grid", k_rule = NULL, k_max = 12),
    list(id = "PCA_BNBIC", factor_method = "PCA", k_mode = "dynamic", k_rule = "bn_bic", k_max = 8)
  )

  expect_true(has_dynamic_rule(specs_mixed, "bn_bic"))
  expect_false(has_dynamic_rule(specs_mixed, "onatski"))

  specs_grid_only <- list(
    list(id = "PCA_grid", factor_method = "PCA", k_mode = "grid", k_rule = NULL, k_max = 12)
  )

  expect_false(has_dynamic_rule(specs_grid_only, "bn_bic"))
})

test_that("filter_specs_by_method works correctly", {
  specs <- list(
    list(id = "PCA_grid", factor_method = "PCA", k_mode = "grid", k_rule = NULL, k_max = 12),
    list(id = "PLS_grid", factor_method = "PLS", k_mode = "grid", k_rule = NULL, k_max = 10),
    list(id = "PCA_BNBIC", factor_method = "PCA", k_mode = "dynamic", k_rule = "bn_bic", k_max = 8)
  )

  pca_specs <- filter_specs_by_method(specs, "PCA")
  expect_length(pca_specs, 2)

  pls_specs <- filter_specs_by_method(specs, "PLS")
  expect_length(pls_specs, 1)
})

test_that("config k_selection_settings has correct defaults", {
  config <- config_us_default()

  expect_true("k_selection_settings" %in% names(config))
  expect_equal(config$k_selection_settings$min_k, 1)
  expect_equal(config$k_selection_settings$bn_bic_sigma_sq, "v_kmax")
  expect_equal(config$k_selection_settings$onatski_r_max, 12)
})

test_that("validate_config validates k_selection_settings", {
  config <- config_us_default()

  # Valid config passes
  expect_silent(validate_config(config))

  # Invalid min_k
  config_bad <- config
  config_bad$k_selection_settings$min_k <- 0
  expect_error(validate_config(config_bad), "min_k must be a positive integer")

  config_bad$k_selection_settings$min_k <- -1
  expect_error(validate_config(config_bad), "min_k must be a positive integer")
})
