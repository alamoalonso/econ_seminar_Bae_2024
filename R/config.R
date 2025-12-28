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

    # File paths
    data_file = "current.csv",
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
    first_forecast_idx = 60,  # Start forecasting at t=60
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
    trace_origins = c(60, 120, 240),  # Time indices to trace for debugging

    # Plotting
    make_plots = TRUE,  # Generate Bae-style figures when saving results

    # Compatibility mode (ensure exact replication)
    compatibility_mode = TRUE
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
  cfg$data_file <- "data/ea-md-qd/2025-11/data_TR2/EAdataM_TR2.xlsx"
  cfg$sample_start <- as.Date("2000-01-01")
  cfg$sample_end <- as.Date("2025-11-01")
  cfg$eval_start <- as.Date("2010-01-01")
  cfg$eval_end <- as.Date("2019-12-01")
  
  cfg$apply_transforms <- FALSE
  
  cfg$first_forecast_idx <- 120

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
