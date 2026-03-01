#' Configuration System
#'
#' Functions to create and manage configuration objects for forecasting workflows.
#'
#' @name config
NULL

#' Create default US FRED-MD configuration
#'
#' @return A list containing all configuration parameters for US dataset
#' @export
#' @examples
#' config <- config_us_default()
config_us_default <- function() {
  list(
    # Dataset identification
    dataset_id = "US_FRED",
    frequency = "monthly",  # Data frequency: "monthly", "quarterly", or "yearly"
                            # Used to convert rolling_window_years to observation count

    # File paths
    data_file = "data/fred-md/current.csv",
    output_dir = "outputs",
    run_id = format(Sys.time(), "%Y%m%d_%H%M%S"),

    # Data parameters
    sample_start = as.Date("1959-03-01"),
    sample_end = as.Date("2019-12-01"),
    eval_start = as.Date("1970-01-01"),
    eval_end = as.Date("2019-12-01"),

    # Series to forecast (NULL = all balanced series)
    # For subset testing, use: c("HWIURATIO","HOUST","CMRMTSPLx","INVEST","AAA","CPIAUCSL","S.P.500")
    series_list = NULL,

    # Forecasting settings
    horizons = c(1, 6, 12, 24),

    # Factor extraction settings
    k_max_pca = 12,
    k_max_pls = 12,
    pca_center = TRUE,
    pca_scale = FALSE,
    pls_center = TRUE,
    pls_scale = FALSE,
    pls_method = "oscorespls",

    # Paper-compliant implementations (Bae 2024 exact specification)
    use_paper_pca = TRUE,  # Use paper_pca_factors() with N^{-1}Λ'Λ = I_r normalization
    use_paper_pls = TRUE,  # Use paper_pls_factors() with Q(F̂_prev)X deflation

    # Model settings
    max_p_ar = 6,        # Max AR lag order for BIC selection
    max_m_diarlag = 3,   # Max factor lag order for DIAR-LAG
    k_max_diarlag = 4,   # Max k for DIAR-LAG (as per Bae 2024)

    # Window settings
    first_forecast_idx = 120,  # Start forecasting at t=120
    rolling_window_years = 10,

    # Models to run
    models = c("AR", "DI", "DIAR", "DIAR-LAG"),
    schemes = c("recursive", "rolling"),
    factor_methods = c("PCA", "PLS"),

    # Preprocessing settings
    apply_transforms = TRUE,
    standardize = TRUE,
    outlier_iqr_multiplier = 10,
    require_balanced_panel = TRUE,

    # Debug and logging
    debug = FALSE,
    suppress_warnings = TRUE,  # Set FALSE to see WARN messages (e.g., Onatski k clamping)
    trace_origins = c(60, 120, 240),  # Time indices to trace for debugging

    # Plotting
    make_plots = TRUE,  # Generate Bae-style figures when saving results

    # Forecast comparison tests
    do_tests = TRUE,  # Compute Diebold-Mariano and Clark-West tests
    test_types = c("DM", "CW"),  # Test types to compute
    test_alpha = c(0.10, 0.05, 0.01),  # Significance levels
    hac_lag_rule = "h-1",  # HAC lag rule: "h-1" means L = max(h-1, 0)

    # Table 5 settings
    table5_factor_method = "PLS",  # Factor method for Table 5
    table5_k = 1,  # Number of factors for Table 5
    category_mapping_file = "data/fred-md/category_mappings.csv",

    # Compatibility mode (ensure exact replication)
    compatibility_mode = TRUE,

    # Forecast persistence for MCS analysis
    save_forecasts = TRUE,  # Persist OOS forecasts at origin level for MCS tests

    # Model Confidence Set (MCS) evaluation
    mcs = list(
      enabled = TRUE,  # Set to TRUE to run MCS tests
      alphas = c(0.10, 0.05, 0.01),  # Significance levels for MCS
      loss = "se",  # Loss function: "se" (squared error) or "ae" (absolute error)
      test_stat = "TR",
      B = 1000,  # Bootstrap replications
      block_length = NULL,  # Block length for bootstrap (NULL = automatic)
      seed = 42,  # Random seed for reproducibility
      include_ar_in_M0 = TRUE,  # If TRUE, always include AR benchmark in MCS comparisons
      # M0_sets: List of candidate method sets to compare
      # Each set is a character vector of method_id patterns or special names:
      # - "k1-PLS", "k1-PCA": PLS/PCA with k=1
      # - "PLS-BNBIC", "PCA-BNBIC": BN-BIC selected k (placeholder)
      # - Regex patterns like "recursive_PCA_DI_k[0-9]+" to match multiple methods
      # Note: AR is automatically included if include_ar_in_M0 = TRUE
      M0_sets = list(
        baseline = c("k1-PLS", "k1-PCA")
      )
    ),

    # Factor extraction specifications
    # If NULL, will be auto-generated from factor_methods for backward compatibility
    # Each spec defines: id, factor_method, k_mode ("grid"/"dynamic"), k_rule, k_max
    # Example with both grid and dynamic methods:
    # factor_specs = list(
    #   list(id = "PCA_grid", factor_method = "PCA", k_mode = "grid", k_rule = NULL, k_max = 12),
    #   list(id = "PLS_grid", factor_method = "PLS", k_mode = "grid", k_rule = NULL, k_max = 12),
    #   list(id = "PCA_BNBIC", factor_method = "PCA", k_mode = "dynamic", k_rule = "bn_bic", k_max = 12),
    #   list(id = "PLS_ON", factor_method = "PLS", k_mode = "dynamic", k_rule = "onatski", k_max = 12)
    # ),
    factor_specs = NULL,  # NULL = auto-generate grid specs from factor_methods

    # k-selection settings (for dynamic rules)
    k_selection_settings = list(
      min_k = 1,                   # Minimum k allowed (enforced: k_hat >= 1)
      bn_bic_sigma_sq = "v_kmax",  # How to compute sigma^2 in BIC3: "v_kmax"
      onatski_r_max = 12,          # r_max for Onatski (auto-adjusted if constraint violated)
      onatski_delta = "default",   # "default" uses max(N^{-2/5}, T^{-2/5})
      fallback_on_error = "k_max"  # What to return if k-selection fails: "k_max" or "min_k"
    )
  )
}

#' Create default Euro Area configuration (stub)
#'
#' @return A list containing all configuration parameters for Euro Area dataset
#' @export
#' @examples
#' config <- config_euro_default()
config_euro_default <- function() {
  # Start with US defaults
  cfg <- config_us_default()

  # Override for Euro Area
  cfg$dataset_id <- "EU_EA_MD_QD"
  cfg$frequency <- "monthly"  # Update to "quarterly" if using quarterly data
  cfg$data_file <- "data/ea-md-qd/2025-11/data_TR2/EAdataM_TR2.xlsx"
  cfg$sample_start <- as.Date("2000-01-01")
  cfg$sample_end <- as.Date("2025-11-01")
  cfg$eval_start <- as.Date("2010-01-01")
  cfg$eval_end <- as.Date("2019-12-01")
  
  cfg$apply_transforms <- FALSE
  
  cfg$first_forecast_idx <- 120
  
  cfg$category_mapping_file <- "data/ea-md-qd/2025-11/category_mappings.csv"

  cfg
}

#' Validate configuration object
#'
#' @param config Configuration list
#' @return Invisibly returns TRUE if valid, stops with error otherwise
#' @export
validate_config <- function(config) {
  required_fields <- c(
    "dataset_id", "data_file", "output_dir",
    "sample_start", "sample_end", "eval_start", "eval_end",
    "horizons", "k_max_pca", "k_max_pls",
    "models", "schemes", "factor_methods"
  )

  missing <- setdiff(required_fields, names(config))
  if (length(missing) > 0) {
    stop("Missing required config fields: ", paste(missing, collapse = ", "))
  }

  # Check date ordering
  if (config$sample_start >= config$sample_end) {
    stop("sample_start must be before sample_end")
  }
  if (config$eval_start < config$sample_start) {
    stop("eval_start cannot be before sample_start")
  }
  if (config$eval_end > config$sample_end) {
    stop("eval_end cannot be after sample_end")
  }

  # Check horizons
  if (!all(config$horizons > 0)) {
    stop("All horizons must be positive integers")
  }

  # Check k_max
  if (config$k_max_pca < 1 || config$k_max_pls < 1) {
    stop("k_max_pca and k_max_pls must be >= 1")
  }

  # Check frequency
  if (!is.null(config$frequency)) {
    valid_frequencies <- c("monthly", "quarterly", "yearly")
    if (!config$frequency %in% valid_frequencies) {
      stop(sprintf("frequency must be one of: %s (got: %s)",
                   paste(valid_frequencies, collapse = ", "),
                   config$frequency))
    }
  }

  # Validate factor_specs if provided
  if (!is.null(config$factor_specs)) {
    validate_factor_specs(config$factor_specs)
  }

  # Validate k_selection_settings if provided
  if (!is.null(config$k_selection_settings)) {
    kss <- config$k_selection_settings
    if (!is.null(kss$min_k) && (kss$min_k < 1 || kss$min_k != as.integer(kss$min_k))) {
      stop("k_selection_settings$min_k must be a positive integer >= 1")
    }
  }

  invisible(TRUE)
}

#' Update configuration with custom parameters
#'
#' @param config Base configuration list
#' @param ... Named arguments to update
#' @return Updated configuration list
#' @export
#' @examples
#' config <- config_us_default()
#' config <- update_config(config, debug = TRUE, horizons = c(1, 12))
update_config <- function(config, ...) {
  updates <- list(...)
  for (nm in names(updates)) {
    config[[nm]] <- updates[[nm]]
  }
  validate_config(config)
  config
}
