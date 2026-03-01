# Factor-Based Forecasting Framework

A modular, package-like R framework for macroeconomic forecasting using factor models (PCA and PLS). Implements the methodology from Bae (2024) with support for multiple forecasting designs: DI (Direct forecast with factors), DIAR (Direct forecast with factors and AR lags), and DIAR-LAG (Direct forecast with factors and their lags).

## Features

- **Modular Architecture**: Clean separation of concerns with dedicated modules for data I/O, preprocessing, factor extraction, forecasting, and evaluation
- **Multiple Factor Methods**: Support for both PCA and PLS factor extraction
- **Multiple Forecasting Designs**: AR, DI, DIAR, DIAR-LAG models
- **Rolling and Recursive Windows**: Both estimation schemes supported
- **Forecast Comparison Tests**: Diebold-Mariano and Clark-West tests with HAC standard errors
- **Model Confidence Set (MCS)**: Hansen-Lunde-Nason (2011) procedure for identifying superior models
- **Forecast Persistence**: Origin-level OOS forecasts saved for post-hoc analysis and MCS
- **Table 5 Replication**: Automatic generation of Bae (2024) Table 5 summaries
- **Performance Optimized**: Unified evaluation computes forecasts only once (~40% faster)
- **Configurable Workflows**: Easy-to-swap datasets and parameters via configuration system
- **Reproducible Results**: Deterministic processing with explicit seed control
- **Comprehensive Logging**: Debug and trace functionality for troubleshooting

## Project Structure

```
.
├── R/                          # Core R modules
│   ├── config.R                # Configuration system (factor_specs, k_selection_settings)
│   ├── data_io.R               # Data loading functions
│   ├── preprocessing.R         # Data transformation and cleaning
│   ├── factors_pca.R           # PCA factor extraction (library-based)
│   ├── factors_pls.R           # PLS factor extraction (library-based)
│   ├── paper_factors_pca.R     # Paper-compliant PCA (Bae 2024 normalization)
│   ├── paper_factors_pls.R     # Paper-compliant PLS
│   ├── factor_specs.R          # Factor spec creation and validation
│   ├── k_selection.R           # BN-BIC and Onatski k-selection rules
│   ├── factor_cache.R          # Per-origin caching for efficient extraction
│   ├── forecasting_designs.R   # Model building/fitting/prediction (AR, DI, DIAR, DIAR-LAG)
│   ├── forecasting_models.R    # High-level forecast orchestration (grid + dynamic)
│   ├── evaluation.R            # RMSE, DM/CW tests, forecast persistence, unified evaluation
│   ├── mcs.R                   # Model Confidence Set (MCS) implementation
│   ├── plots_rmse.R            # Bae-style RMSE plotting functions
│   ├── workflow.R              # Top-level workflow orchestration
│   └── utils_logging.R         # Logging utilities
├── scripts/
│   ├── run_workflow.R          # Main entry point script
│   ├── example_usage.R         # Example usage patterns
│   └── validate_optimization.R # Validation script for performance optimization
├── inst/extdata/               # Optional: example config files, sample data
├── tests/testthat/             # Unit and smoke tests
│   ├── test_smoke.R            # Integration tests
│   ├── test_k_selection.R      # k-selection unit tests
│   ├── test_factor_specs.R     # Factor specs unit tests
│   └── ...
├── outputs/                    # Results output directory (created automatically)
├── DESCRIPTION                 # Package metadata
├── NAMESPACE                   # Package namespace
└── README.md                   # This file
```

## Installation

This project is structured as an R package but does not need to be formally installed. Simply ensure all dependencies are available:

```r
# Install required packages
install.packages(c(
  "fbi",        # For FRED-MD data handling
  "pls",        # For PLS factor extraction
  "tidyverse",  # Data manipulation and visualization
  "dplyr",
  "tibble",
  "readr",
  "ggplot2",
  "patchwork",  # For combining plots with independent scales
  "purrr",
  "lubridate",
  "readxl",
  "janitor"
))

# Optional: For efficient forecast persistence (parquet format)
install.packages("arrow")

# Recommended: For Model Confidence Set evaluation
install.packages("MCS")
```

## Quick Start

### Running the US Workflow

The simplest way to run the complete workflow:

```bash
# From the project root directory
Rscript scripts/run_workflow.R US
```

This will:
1. Load FRED-MD data from `current.csv`
2. Preprocess and standardize the data
3. Extract PCA and PLS factors
4. Run forecasts for all series and horizons
5. Compute RMSE and forecast comparison tests (DM/CW) - **optimized to compute forecasts only once**
6. Save origin-level OOS forecasts for MCS analysis (default: enabled)
7. Compute MCS evaluation (if `config$mcs$enabled = TRUE`)
8. Generate Table 5 summaries (if category mapping provided)
9. Save results to `outputs/[timestamp]/`

### Running from R Console

```r
# Set working directory to project root
setwd("path/to/project")

# Source all modules
source_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (f in source_files) source(f)

# Load required libraries
library(fbi)
library(pls)
library(tidyverse)
library(lubridate)

# Create configuration
config <- config_us_default()

# Optional: Provide category mapping for Table 5 generation
# config$category_mapping_file <- "data/fred-md/category_mappings.csv"

# Optional: Enable MCS evaluation
# config$mcs$enabled <- TRUE

# Run workflow
results <- run_workflow(config)

# Compute evaluation (RMSE + tests + forecasts) - OPTIMIZED
# This computes forecasts once and derives both RMSE and tests
evaluation <- compute_evaluation(results, config)
rmse_results <- evaluation$rmse_results
tests_results <- evaluation$tests_results
forecasts <- evaluation$forecasts  # Origin-level OOS forecasts

# Optional: Compute MCS if enabled
mcs_results <- NULL
if (isTRUE(config$mcs$enabled) && !is.null(forecasts)) {
  mcs_results <- compute_mcs_evaluation(forecasts, config)
}

# Save results (includes forecasts, MCS results, tests, Table 5)
save_results(results, rmse_results, config, tests_results,
             forecasts = forecasts, mcs_results = mcs_results)

# Print summary
print_summary(results, rmse_results, config, tests_results)
```

## Configuration System

The framework uses a configuration-based approach to manage all parameters. Two default configurations are provided:

### US FRED-MD Configuration

```r
config <- config_us_default()
```

Key parameters:
- `dataset_id`: "US_FRED"
- `data_file`: "current.csv"
- `sample_start`: 1959-03-01
- `sample_end`: 2019-12-01
- `horizons`: c(1, 6, 12, 24)
- `k_max_pca`: 12
- `k_max_pls`: 12

### Euro Area Configuration (Stub)

```r
config <- config_euro_default()
```

This is a placeholder configuration. Update the following:
- `data_file`: Path to your Euro Area CSV file
- `sample_start/end`: Appropriate date range
- Other parameters as needed

### Customizing Configuration

```r
# Start with default
config <- config_us_default()

# Modify specific parameters
config <- update_config(
  config,
  debug = TRUE,
  horizons = c(1, 12),
  series_list = c("INDPRO", "CPIAUCSL", "S.P.500")
)

# Run with custom config
results <- run_workflow(config)
```

### Key Configuration Options

#### Data and Sample Settings
- `dataset_id`: Dataset identifier ("US_FRED" or "EURO_AREA")
- `data_file`: Path to input CSV file
- `sample_start/end`: Sample period boundaries
- `eval_start/end`: Evaluation period boundaries
- `series_list`: Specific series to forecast (NULL = all balanced series)

#### Model and Forecasting Settings
- `horizons`: Forecast horizons (e.g., c(1, 6, 12, 24))
- `k_max_pca/pls`: Maximum number of factors to extract
- `models`: Models to run (c("AR", "DI", "DIAR", "DIAR-LAG"))
- `schemes`: Estimation schemes (c("recursive", "rolling"))
- `factor_methods`: Factor extraction methods (c("PCA", "PLS"))

#### Forecast Comparison Tests
- `do_tests`: Enable/disable test computation (default: TRUE)
- `test_types`: Which tests to compute (default: c("DM", "CW"))
- `test_alpha`: Significance levels (default: c(0.10, 0.05, 0.01))
- `hac_lag_rule`: HAC lag specification (default: "h-1" meaning L = max(h-1, 0))

#### Table 5 Settings
- `table5_factor_method`: Factor method for Table 5 (default: "PLS")
- `table5_k`: Number of factors for Table 5 (default: 1)
- `category_mapping_file`: Path to CSV with series,category columns (required for Table 5)

#### Forecast Persistence
- `save_forecasts`: Persist OOS forecasts at origin level for MCS analysis (default: TRUE)

#### Model Confidence Set (MCS)
- `mcs$enabled`: Enable MCS evaluation (default: FALSE)
- `mcs$alphas`: Significance levels (default: c(0.10, 0.05))
- `mcs$loss`: Loss function - "se" (squared error) or "ae" (absolute error)
- `mcs$test_stat`: Test statistic - "TR" (range) or "Tmax" (max)
- `mcs$B`: Bootstrap replications (default: 1000)
- `mcs$seed`: Random seed for reproducibility
- `mcs$M0_sets`: Named list of candidate method sets to compare

#### Debugging and Output
- `debug`: Enable debug logging (TRUE/FALSE)
- `trace_origins`: Time indices to trace for detailed logging
- `make_plots`: Generate Bae-style figures (default: TRUE)

## Forecast Comparison Tests

The framework implements one-sided Diebold-Mariano (DM) and Clark-West (CW) tests to compare competing models against the AR benchmark.

### Test Specifications

**Diebold-Mariano (DM) Test:**
- Loss differential: `d_t = e0_t^2 - e1_t^2` where e0 = AR error, e1 = competing model error
- Null hypothesis: H0: E[d_t] ≤ 0
- Alternative hypothesis: H1: E[d_t] > 0 (competing model is more accurate)
- Standard errors: Newey-West HAC with lag L = max(h-1, 0) and Bartlett kernel

**Clark-West (CW) Test:**
- Adjusted differential: `d_t = e0_t^2 - (e1_t^2 - delta_f_t^2)` where delta_f = f0 - f1
- Accounts for nesting bias in nested models
- Same hypotheses and HAC specification as DM test

### Configuration

```r
config <- config_us_default()

# Enable/disable tests
config$do_tests <- TRUE  # default: TRUE

# Which tests to compute
config$test_types <- c("DM", "CW")  # default: both

# Significance levels
config$test_alpha <- c(0.10, 0.05, 0.01)  # default

# HAC lag rule
config$hac_lag_rule <- "h-1"  # L = max(h-1, 0)
```

### Interpreting Test Results

```r
# Load test results
tests <- read.csv("outputs/[run_id]/tests_results.csv")

# Filter significant results at 5% level
significant <- tests %>%
  filter(reject_005 == TRUE) %>%
  select(series_id, h, scheme, model_class, k, test_type, stat, p_value_one_sided)

# Count rejections by model
rejection_summary <- tests %>%
  group_by(model_class, test_type, scheme) %>%
  summarise(
    total_tests = n(),
    rejections_10pct = sum(reject_010, na.rm = TRUE),
    rejections_5pct = sum(reject_005, na.rm = TRUE),
    rejections_1pct = sum(reject_001, na.rm = TRUE)
  )
```

## Category Mapping for Table 5

To generate Bae (2024) Table 5 summaries, you must provide a CSV file mapping each series to its category.

### Creating the Category Mapping File

Create a CSV file (e.g., `data/category_mapping.csv`) with two columns:

```csv
series,category
RPI,Output and Income
W875RX1,Output and Income
INDPRO,Output and Income
IPFPNSS,Output and Income
IPFINAL,Output and Income
IPCONGD,Output and Income
IPDCONGD,Output and Income
CUMFNS,Employment and Hours
HWI,Employment and Hours
HWIURATIO,Employment and Hours
CLF16OV,Employment and Hours
CE16OV,Employment and Hours
UNRATE,Employment and Hours
UEMPMEAN,Employment and Hours
...
```

**Category examples from Bae (2024):**
- Output and Income
- Employment and Hours
- Consumption and Orders
- Inventories and Sales
- Prices
- Interest Rates and Spreads
- Money and Credit
- Stock Market

### Using Category Mapping

```r
# Set category mapping file in config
config <- config_us_default()
config$category_mapping_file <- "data/category_mapping.csv"

# Run workflow
results <- run_workflow(config)
evaluation <- compute_evaluation(results, config)

# Save results - Table 5 will be generated automatically
save_results(results, evaluation$rmse_results, config, evaluation$tests_results)
```

The workflow will automatically:
1. Load the category mapping
2. Join categories to test results
3. Compute rejection frequencies by category
4. Generate `table5_dm.csv` and `table5_cw.csv`

**Important**: All series in your dataset must be present in the category mapping file, or Table 5 generation will fail with a clear error message.

## Model Confidence Set (MCS)

The framework integrates with the [MCS package](https://cran.r-project.org/package=MCS) (Catania & Bernardi) to implement the Hansen-Lunde-Nason (2011) Model Confidence Set procedure for identifying superior forecasting methods.

**Important**: Install the MCS package for correct results: `install.packages("MCS")`

### What is MCS?

The Model Confidence Set is a sequential testing procedure that identifies a set of models that are statistically indistinguishable from the best model at a given significance level. Unlike pairwise tests (DM/CW), MCS provides:

- **Multiple comparison correction**: Controls for the family-wise error rate
- **Superior set identification**: Returns the set of models that cannot be rejected as inferior
- **Bootstrap inference**: Uses block bootstrap with AR-based block length selection

### Enabling MCS

```r
config <- config_us_default()

# Enable MCS evaluation
config$mcs$enabled <- TRUE

# Configure MCS parameters
config$mcs$alphas <- c(0.10, 0.05)  # Test at multiple significance levels
config$mcs$loss <- "se"             # "se" (squared error) or "ae" (absolute error)
config$mcs$test_stat <- "TR"        # "TR" (range statistic) or "Tmax" (max statistic)
config$mcs$B <- 1000                # Bootstrap replications
config$mcs$seed <- 42               # For reproducibility

# Include AR benchmark in all MCS comparisons (default: TRUE)
# When TRUE, AR is automatically added to each M0 comparison
config$mcs$include_ar_in_M0 <- TRUE

# Define candidate method sets to compare
# Note: AR is automatically included if include_ar_in_M0 = TRUE
config$mcs$M0_sets <- list(
  # Compare k=1 factor methods (AR included automatically)
  baseline = c("k1-PLS", "k1-PCA"),

  # Compare all PCA-based methods (AR included automatically)
  all_pca = c("k1-PCA", "k2-PCA", "k3-PCA", "k4-PCA")
)

# Run workflow
results <- run_workflow(config)
evaluation <- compute_evaluation(results, config)

# Compute MCS
mcs_results <- compute_mcs_evaluation(evaluation$forecasts, config)

# Save all results
save_results(results, evaluation$rmse_results, config, evaluation$tests_results,
             forecasts = evaluation$forecasts, mcs_results = mcs_results)
```

**Note on slicing:** MCS tests are run separately for each forecast equation (model_class: DI, DIAR, DIAR-LAG). Within each slice, the specified factor methods are compared, and if `include_ar_in_M0 = TRUE`, the AR benchmark is included to test whether factor-based methods are statistically distinguishable from the benchmark.

### Method Specification Syntax

The `M0_sets` parameter accepts various method specifications:

| Syntax | Description | Example Match |
|--------|-------------|---------------|
| `"AR"` | AR benchmark | `recursive_AR`, `rolling_AR` |
| `"k1-PCA"` | PCA with k=1 factors | `recursive_PCA_DI_k1`, `rolling_PCA_DIAR_k1` |
| `"k3-PLS"` | PLS with k=3 factors | `recursive_PLS_DI_k3`, `rolling_PLS_DIAR_k3` |
| Regex pattern | Custom matching | `"recursive_PCA_.*_k[1-4]"` |

### Interpreting MCS Results

```r
# Load MCS results
mcs <- arrow::read_parquet("outputs/[run_id]/mcs_results.parquet")
# Or: mcs <- read.csv("outputs/[run_id]/mcs_results.csv")

# View superior sets at alpha = 0.10
mcs %>%
  filter(alpha == 0.10, status == "ok") %>%
  select(series_id, h, scheme, M0_set_id, superior_set, n_methods_mcs)

# Count how often each method is in the MCS
library(tidyr)
mcs %>%
  filter(status == "ok") %>%
  separate_rows(superior_set, sep = ";") %>%
  count(superior_set, sort = TRUE)
```

### MCS Output Structure

The `mcs_results.parquet` (or `.csv`) contains:

| Column | Description |
|--------|-------------|
| `run_id`, `series_id`, `h`, `scheme`, `model_class` | Slice identifiers |
| `alpha` | Significance level |
| `M0_set_id` | Name of the candidate method set |
| `loss` | Loss function used |
| `test_stat` | Test statistic type |
| `B` | Bootstrap replications |
| `n_methods_input` | Number of methods in candidate set |
| `n_methods_mcs` | Number of methods in the MCS |
| `included_methods` | Semicolon-separated input methods |
| `superior_set` | Semicolon-separated MCS members |
| `pvalues` | Method=p-value pairs |
| `status` | "ok", "skipped", or "error" |

### Building Loss Matrices Manually

For custom MCS analysis, you can build loss matrices from saved forecasts:

```r
# Build loss matrix for a specific series and horizon
L_result <- build_mcs_loss_matrix(
  forecasts_path = "outputs/[run_id]/forecasts",
  series_id = "INDPRO",
  h = 1,
  loss_fn = function(y, yhat) (y - yhat)^2
)

# L_result contains:
# - L: Matrix (T_oos x M) of losses
# - methods: Vector of method names
# - origin_dates: Dates for each row
# - n_dropped: Rows dropped due to NA (MCS requires balanced panel)

# Run MCS manually
mcs_result <- run_mcs(L_result$L, alpha = 0.10, B = 1000, stat_type = "TR")
print(mcs_result$superior_set)
print(mcs_result$pvalues)
```

### References

**Hansen, P. R., Lunde, A., & Nason, J. M. (2011)**. The Model Confidence Set. *Econometrica*, 79(2), 453-497.

## Performance Optimization

The framework uses a **unified evaluation** approach that computes forecasts only once, then derives both RMSE and test statistics from the same forecast objects.

### Optimized Approach (Recommended)

```r
# Single unified call - computes forecasts ONCE
evaluation <- compute_evaluation(results, config)
rmse_results <- evaluation$rmse_results
tests_results <- evaluation$tests_results
forecasts <- evaluation$forecasts  # Origin-level OOS forecasts (if save_forecasts=TRUE)
```

**Performance gain**: ~40% faster than computing RMSE and tests separately

**Bonus**: When `save_forecasts = TRUE` (default), origin-level forecasts are collected for MCS analysis without recomputing.

### Legacy Approach (Backward Compatible)

```r
# Separate calls - computes forecasts TWICE (less efficient)
rmse_results <- compute_rmse(results, config)
tests_results <- compute_tests(results, config)
```

Use this only if you need RMSE or tests alone, or for backward compatibility.

### Validation

To verify the optimization produces identical results:

```bash
Rscript scripts/validate_optimization.R
```

This will:
- Run both approaches on a small dataset
- Compare results (should be numerically identical)
- Measure performance improvement
- Extrapolate time savings to full dataset

See `OPTIMIZATION_NOTES.md` for detailed performance analysis.

## Swapping Datasets

### Method 1: Create New Configuration

```r
# Define custom configuration
config_custom <- config_us_default()
config_custom$dataset_id <- "MY_DATASET"
config_custom$data_file <- "my_data.csv"
config_custom$sample_start <- as.Date("2000-01-01")
config_custom$sample_end <- as.Date("2023-12-01")

# Ensure CSV has 'date' column + predictor columns
results <- run_workflow(config_custom)
```

### Method 2: Extend Data Loader

Add a new loader function in `R/data_io.R`:

```r
load_my_dataset <- function(config) {
  # Load your data
  df <- read.csv(config$data_file)
  df$date <- as.Date(df$date)

  # Apply transformations as needed
  # ...

  list(
    raw_data = df,
    dates = df$date,
    dataset_id = config$dataset_id
  )
}
```

Update `load_dataset()` to dispatch to your loader:

```r
load_dataset <- function(config) {
  if (config$dataset_id == "MY_DATASET") {
    return(load_my_dataset(config))
  }
  # ... existing code ...
}
```

## Output Structure

Results are saved to `outputs/[run_id]/`:

```
outputs/20231214_153045/
├── rmse_results.csv                    # Full RMSE table
├── tests_results.csv                   # Forecast comparison test results (DM/CW)
├── table5_dm.csv                       # Table 5 summary (Diebold-Mariano)
├── table5_cw.csv                       # Table 5 summary (Clark-West)
├── forecasts/                          # OOS forecasts for MCS analysis (partitioned parquet)
│   └── series_id=INDPRO/h=1/...        # Partitioned by series and horizon
├── config.rds                          # Configuration used for this run
├── summary.txt                         # Text summary of the run
├── fig1_mean_rmse_k_pca.png            # PCA RMSE plot (mean across all series)
├── fig2_mean_rmse_k_pls.png            # PLS RMSE plot (mean across all series)
├── fig1_mean_rmse_k_pca_INDPRO.png     # PCA RMSE plot for INDPRO series
├── fig1_mean_rmse_k_pca_CPIAUCSL.png   # PCA RMSE plot for CPIAUCSL series
├── fig2_mean_rmse_k_pls_INDPRO.png     # PLS RMSE plot for INDPRO series
└── ...                                 # Additional series-specific plots
```

**Note**: If the `arrow` package is not installed, forecasts are saved as `forecasts_long.csv` instead of partitioned parquet.

### RMSE Results Table

The `rmse_results.csv` contains:
- `series`: Series name
- `h`: Forecast horizon
- `scheme`: "recursive" or "rolling"
- `factor_method`: "PCA" or "PLS"
- `model`: "AR", "DI", "DIAR", or "DIAR-LAG"
- `k`: Number of factors used
- `mse`: Mean squared error
- `mse_ar`: MSE of AR benchmark
- `rmse_rel`: RMSE relative to AR (sqrt(mse / mse_ar))

### Test Results Table

The `tests_results.csv` contains one row per test with:
- `series_id`: Series name
- `scheme`: "recursive" or "rolling"
- `factor_method`: "PCA" or "PLS"
- `model_class`: "DI", "DIAR", or "DIAR-LAG"
- `k`: Number of factors used
- `h`: Forecast horizon
- `test_type`: "DM" or "CW"
- `stat`: Test statistic
- `p_value_one_sided`: One-sided p-value (H1: competing model is more accurate)
- `n_oos`: Number of out-of-sample observations
- `reject_010`, `reject_005`, `reject_001`: Rejection indicators at α = 0.10, 0.05, 0.01

### Table 5 Summaries

The `table5_dm.csv` and `table5_cw.csv` files replicate Bae (2024) Table 5, containing:
- `category`: Series category (e.g., "Labor Market", "Prices", "Overall")
- `scheme`: "recursive" or "rolling"
- `alpha`: Significance level (0.10, 0.05, 0.01)
- `test_type`: "DM" or "CW"
- `frequency`: Number of rejections in category
- `total`: Total number of tests in category
- `percentage`: Percentage of rejections (100 × frequency / total)

**Note**: Table 5 generation requires a category mapping file (see Category Mapping section below).

### OOS Forecasts for MCS Analysis

The `forecasts/` directory (or `forecasts_long.csv`) contains origin-level out-of-sample forecasts for Model Confidence Set (MCS) analysis:

- `run_id`: Run identifier
- `series_id`: Series name
- `h`: Forecast horizon
- `scheme`: "recursive" or "rolling"
- `factor_method`: "PCA", "PLS", or NA (for AR model)
- `model_class`: "AR", "DI", "DIAR", or "DIAR-LAG"
- `k`: Number of factors (NA for AR)
- `origin_index`: Time index at forecast origin
- `origin_date`: Date at forecast origin
- `target_index`: Time index of target realization (origin_index + h)
- `target_date`: Date of target realization
- `y_true`: Realized value
- `y_hat`: Forecast value
- `method_id`: Unique method identifier for MCS (e.g., "recursive_PCA_DI_k3")

### MCS Results Table

When MCS evaluation is enabled, `mcs_results.parquet` (or `.csv`) contains:

- `run_id`, `series_id`, `h`, `scheme`, `model_class`: Slice identifiers
- `alpha`: Significance level
- `M0_set_id`: Name of the candidate method set
- `loss`: Loss function used ("se" or "ae")
- `test_stat`: Test statistic ("TR" or "Tmax")
- `B`: Number of bootstrap replications
- `n_methods_input`: Number of methods in the candidate set
- `n_methods_mcs`: Number of methods in the Model Confidence Set
- `included_methods`: Semicolon-separated list of input methods
- `superior_set`: Semicolon-separated list of methods in the MCS
- `pvalues`: Semicolon-separated list of method=p-value pairs
- `status`: "ok", "skipped", or "error"
- `message`: Status message

See the [Model Confidence Set (MCS)](#model-confidence-set-mcs) section for configuration details.

### Forecast Persistence Configuration

```r
config <- config_us_default()
config$save_forecasts <- TRUE  # Default: persist forecasts for MCS
config$save_forecasts <- FALSE # Disable to save disk space
```

## Plotting

The framework includes plotting functions to replicate the publication-quality figures from Bae (2024).

### Automatic Plot Generation

Enable automatic plot generation during the workflow:

```r
config <- config_us_default()
config$make_plots <- TRUE

results <- run_workflow(config)
rmse_results <- compute_rmse(results, config)
save_results(results, rmse_results, config)  # Plots saved automatically
```

This will generate:
- **Overall mean plots** (averaged across all series):
  - `fig1_mean_rmse_k_pca.png`: Mean RMSE of k-PCA across schemes and models
  - `fig2_mean_rmse_k_pls.png`: Mean RMSE of k-PLS across schemes and models
- **Individual series plots** (one plot per series):
  - `fig1_mean_rmse_k_pca_SERIESNAME.png`: k-PCA results for each specific series
  - `fig2_mean_rmse_k_pls_SERIESNAME.png`: k-PLS results for each specific series
  - Each plot has a dynamic title indicating the series name (e.g., "Mean RMSE of k-PCA for series INDPRO")
  - The set of series plotted is controlled by `config$series_list` (defaults to all series if NULL)

### Manual Plot Generation

Generate plots manually from existing RMSE results:

```r
# Load existing results and config
rmse_results <- read.csv("outputs/[run_id]/rmse_results.csv")
config <- readRDS("outputs/[run_id]/config.rds")

# Generate Bae-style figures (generates both overall and per-series plots)
output_dir <- "outputs/[run_id]"
plot_bae_fig1_kpca(rmse_results, metric = "rmse_rel", save_dir = output_dir, config = config)
plot_bae_fig2_kpls(rmse_results, metric = "rmse_rel", save_dir = output_dir, config = config)
```

**Note**: These functions will generate:
1. Overall mean plots averaged across all series
2. Individual plots for each series specified in `config$series_list` (or all unique series in `rmse_results` if `config$series_list` is NULL)

### Generate Plots for Already-Completed Runs

If you have runs that were saved without plots (e.g., before `make_plots` was available), you can easily generate plots for them:

```r
# Generate plots for a specific completed run
generate_plots_for_run("outputs/20231214_153045")

# Generate plots for the most recent run
runs <- list.dirs("outputs", recursive = FALSE)
latest_run <- runs[length(runs)]
generate_plots_for_run(latest_run)

# Generate plots for all completed runs at once
all_runs <- list.dirs("outputs", recursive = FALSE)
generate_plots_for_multiple_runs(all_runs)

# Generate plots with absolute RMSE instead of relative
generate_plots_for_run("outputs/20231214_153045", metric = "rmse")
```

This is particularly useful when:
- You ran the workflow before plotting functionality was added
- You want to regenerate plots with a different metric (absolute vs relative RMSE)
- You need to batch-generate plots for multiple historical runs

### Custom Plotting

Create custom plots for specific factor methods or metrics:

```r
# Plot PCA results with absolute RMSE instead of relative
p <- plot_mean_rmse_kfactor(
  rmse_tbl = rmse_results,
  factor_method = "PCA",
  metric = "rmse",  # absolute RMSE
  schemes = c("recursive", "rolling"),
  models = c("DI", "DIAR", "DIAR-LAG"),
  horizons = c(1, 6, 12, 24),
  save_path = "my_custom_plot.png",
  width = 10,
  height = 9,
  dpi = 300
)

# Or return plot object without saving
p <- plot_mean_rmse_kfactor(
  rmse_tbl = rmse_results,
  factor_method = "PLS",
  metric = "rmse_rel",
  save_path = NULL  # Don't save, just return plot object
)

# Customize further with ggplot2
library(ggplot2)
p + labs(title = "My Custom Title") + theme_bw()
```

### Data Summarization

Prepare summarized data for plotting or analysis:

```r
# Compute mean RMSE by k across series
summary_data <- summarise_mean_rmse_by_k(
  rmse_tbl = rmse_results,
  metric = "rmse_rel",
  models = c("DI", "DIAR", "DIAR-LAG"),
  horizons = c(1, 6, 12, 24),
  k_max_diarlag = 4
)

# View or export summary
print(summary_data)
write.csv(summary_data, "rmse_summary_by_k.csv", row.names = FALSE)
```

### Plot Structure

The generated figures replicate Bae (2024) Figures 1 and 2:
- **2×3 panel layout**: 2 rows (recursive, rolling) × 3 columns (DI, DIAR, DIAR-LAG)
- **X-axis**: Number of factors k (1-12 for DI/DIAR, 1-4 for DIAR-LAG)
- **Y-axis**: Mean RMSE across all series (relative to AR benchmark by default)
  - **Independent scales per model**: Each model column (DI, DIAR, DIAR-LAG) has its own y-axis scale
  - Implementation: Uses `facet_wrap(~ model, scales = "free_y")` for each scheme, then combines vertically with `patchwork`
  - This prevents DI's large k=1 values from compressing the DIAR/DIAR-LAG panels
- **Lines**: Different forecast horizons (h = 1, 6, 12, 24) shown with distinct line styles and markers
- **Legend**: Single shared legend at bottom combining all horizon information
- **Style**: Grayscale theme with minimal design for publication quality

## Advanced Usage

### Running a Subset of Series

```r
config <- config_us_default()
config$series_list <- c("INDPRO", "CPIAUCSL", "S.P.500")
results <- run_workflow(config)
```

Or via environment variable:

```bash
SERIES_SUBSET="INDPRO,CPIAUCSL,S.P.500" Rscript scripts/run_workflow.R
```

### Enabling Debug Mode

```r
config <- config_us_default()
config$debug <- TRUE
config$trace_origins <- c(60, 120, 240)  # Trace these time indices
results <- run_workflow(config)
```

Or via environment variable:

```bash
DEBUG=true Rscript scripts/run_workflow.R
```

### Running Individual Components

```r
# Load and preprocess data only
config <- config_us_default()
data <- load_dataset(config)
dataset <- preprocess_dataset(data, config)

# Extract factors for a specific origin
dataset$targets_list <- construct_targets(dataset$panel_final, config$horizons)
t_origin <- dataset$dates[120]
fcts <- extract_factors_at_origin(
  panel_final = dataset$panel_final,
  targets_list = dataset$targets_list,
  target_name = "INDPRO",
  h = 1,
  t_origin = t_origin,
  k_max_pca = 12,
  k_max_pls = 12,
  config = config
)

# Forecast a single series
forecasts <- run_forecasts_for_series(
  series_name = "INDPRO",
  h = 1,
  dates = dataset$dates,
  panel_final = dataset$panel_final,
  panel_std1 = dataset$panel_std1,
  targets_list = dataset$targets_list,
  factor_method = "PCA",
  k_max_pca = 12,
  k_max_pls = 12,
  config = config
)
```

## Model Descriptions

### AR (Autoregressive)
- Benchmark model using only lagged values of the target series
- Lag order selected via BIC (max 6 lags)

### DI (Direct Forecast with Factors)
- Direct h-step forecast using k factors
- Regression: y_{t+h} ~ const + F_t

### DIAR (Direct Forecast with Factors and AR Lags)
- Augments DI with AR lags of the target series
- Regression: y_{t+h} ~ const + F_t + y_{t-1} + ... + y_{t-p}
- AR lag order p selected via BIC (max 6)

### DIAR-LAG (Direct Forecast with Factors and Their Lags)
- Uses current and lagged factors (k ≤ 4 per Bae 2024)
- Regression: y_{t+h} ~ const + F_t + F_{t-1} + ... + F_{t-m} + y_{t-1} + ... + y_{t-p}
- Both p and m selected via BIC (max p=6, max m=3)

## Factor Number Selection (k-Selection Rules)

The framework supports both **fixed grid** (evaluate k=1..k_max) and **dynamic** (select k_hat per origin) factor number selection in the same workflow run.

### Available Decision Rules

| Rule | Method | Description |
|------|--------|-------------|
| **Fixed Grid** | PCA, PLS | Evaluate all k from 1 to k_max |
| **BN-BIC** | PCA | Bai-Ng (2002) BIC3 criterion selects k_hat at each origin |
| **Onatski** | PLS | Onatski (2010) eigenvalue edge detection selects k_hat at each origin |

### Configuration

#### Using Factor Specs (Recommended)

Define explicit factor specifications to run both grid and dynamic methods in one workflow:

```r
config <- config_us_default()

# Define mixed grid + dynamic specs
config$factor_specs <- list(
  # Grid methods: evaluate k = 1..k_max
  create_factor_spec("PCA_grid", "PCA", "grid", NULL, 12),
  create_factor_spec("PLS_grid", "PLS", "grid", NULL, 12),

  # Dynamic methods: select k_hat(t) at each origin
  create_factor_spec("PCA_BNBIC", "PCA", "dynamic", "bn_bic", 12),
  create_factor_spec("PLS_ON", "PLS", "dynamic", "onatski", 12)
)

# k-selection settings
config$k_selection_settings <- list(
  min_k = 1,                   # Enforce k_hat >= 1 (no zero factors)
  bn_bic_sigma_sq = "v_kmax",  # sigma^2 = V(k_max) for BIC3
  onatski_r_max = 12           # r_max for Onatski edge detection
)

# Run workflow with mixed specs
results <- run_workflow(config)
evaluation <- compute_evaluation_with_specs(results, config)
```

#### Backward Compatibility

If `factor_specs` is NULL, the framework auto-generates grid specs from `factor_methods`:

```r
config <- config_us_default()
config$factor_methods <- c("PCA", "PLS")  # Auto-expands to PCA_grid, PLS_grid
```

### Factor Spec Structure

Each factor spec has the following fields:

| Field | Values | Description |
|-------|--------|-------------|
| `id` | String | Unique identifier (e.g., "PCA_grid", "PLS_ON") |
| `factor_method` | "PCA", "PLS", "1-PLS" | Factor extraction method |
| `k_mode` | "grid", "dynamic" | Grid evaluates k=1..k_max; dynamic computes k_hat |
| `k_rule` | NULL, "bn_bic", "onatski" | Decision rule (only for dynamic mode) |
| `k_max` | Integer | Maximum k (upper bound for both modes) |

### BN-BIC (Bai-Ng BIC3) for PCA

At each forecast origin t with training data X (T_t × N):

```
For k in {1, ..., k_max}:
    V(k) = (1/(N*T_t)) * ||X - X_k||²   # Reconstruction residual variance

σ² = V(k_max)   # Default convention

BIC3(k) = ln(V(k)) + k * σ² * ((N + T_t - k)/(N*T_t)) * ln(N*T_t)

k_hat = argmin BIC3(k), subject to k_hat >= min_k
```

**Implementation**: Uses single SVD decomposition for efficiency (V(k) computed from cumulative eigenvalues).

### Onatski (2010) for PLS

At each forecast origin t with training data X (T_t × N):

```
1. Compute eigenvalues λ_1 ≥ ... ≥ λ_N of (1/T_t) X'X
2. Set r_max = min(config$onatski_r_max, floor((min(N,T_t)-1)/2))
3. Compute w = 2^(2/3) / (2^(2/3) - 1) ≈ 2.7321
4. Compute u_hat = w * λ_{r_max+1} + (1-w) * λ_{2*r_max+1}
5. Compute δ = max(N^{-2/5}, T_t^{-2/5})
6. Select k_hat = #{i : λ_i > (1+δ) * u_hat}, subject to k_hat >= min_k
```

**Important**: Onatski eigenvalues are computed on PLS-preprocessed X (same center/scale as PLS extraction) for consistency.

### Output Schema

The enhanced forecast output includes k-selection metadata:

| Column | Description |
|--------|-------------|
| `factor_spec_id` | Spec identifier (e.g., "PCA_grid", "PLS_ON") |
| `k_mode` | "grid" or "dynamic" |
| `k_selection_rule` | "fixed" for grid; "bn_bic" or "onatski" for dynamic |
| `k` | Actual k used (1..k_max for grid; k_hat for dynamic) |
| `training_window_start` | First index in training set |
| `training_window_end` | Last index in training set |

### Method ID Convention

| Type | Format | Example |
|------|--------|---------|
| Grid | `{scheme}_{spec_id}_{model}_k{n}` | `recursive_PCA_grid_DI_k3` |
| Dynamic | `{scheme}_{spec_id}_{model}` | `recursive_PCA_BNBIC_DI` |
| AR | `{scheme}_AR` | `recursive_AR` |

**Note**: Dynamic methods have no `_k{n}` suffix because k varies per origin.

### MCS Method Resolution

The `resolve_method_specs()` function supports both old and new formats:

```r
# Match dynamic BN-BIC methods
resolve_method_specs("PCA_BNBIC", available_methods, scheme = "recursive")
# Returns: "recursive_PCA_BNBIC_DI", "recursive_PCA_BNBIC_DIAR", ...

# Match dynamic Onatski methods
resolve_method_specs("PLS_ON", available_methods, scheme = "rolling")
# Returns: "rolling_PLS_ON_DI", "rolling_PLS_ON_DIAR", ...

# Match grid methods
resolve_method_specs("PCA_grid", available_methods, model_class = "DI")
# Returns: "recursive_PCA_grid_DI_k1", ..., "recursive_PCA_grid_DI_k12"
```

### Files and Functions

| File | Key Functions |
|------|---------------|
| `R/factor_specs.R` | `create_factor_spec()`, `get_factor_specs()`, `validate_factor_specs()` |
| `R/k_selection.R` | `compute_bn_bic_k()`, `compute_onatski_k()`, `select_k_dynamic()` |
| `R/forecasting_models.R` | `run_forecasts_for_spec()`, `run_forecasts_for_all_specs()` |
| `R/evaluation.R` | `compute_evaluation_with_specs()`, `extract_forecasts_from_spec_res()` |

## Testing

### Smoke Test

Run a quick smoke test on a small sample:

```r
config <- config_us_default()
config$series_list <- c("INDPRO")
config$horizons <- c(1)
config$factor_methods <- c("PCA")
config$do_tests <- TRUE
config$save_forecasts <- TRUE  # Default, but explicit for clarity

results <- run_workflow(config)
evaluation <- compute_evaluation(results, config)

# Verify output structure
stopifnot(nrow(evaluation$rmse_results) > 0)
stopifnot(all(c("series", "h", "scheme", "model", "mse") %in% names(evaluation$rmse_results)))

# Verify test results if enabled
if (!is.null(evaluation$tests_results)) {
  stopifnot(nrow(evaluation$tests_results) > 0)
  stopifnot(all(c("series_id", "test_type", "stat", "p_value_one_sided") %in% names(evaluation$tests_results)))
}

# Verify forecasts are collected (for MCS analysis)
if (!is.null(evaluation$forecasts)) {
  stopifnot(nrow(evaluation$forecasts) > 0)
  stopifnot(all(c("series_id", "h", "method_id", "y_true", "y_hat") %in% names(evaluation$forecasts)))
}
```

### Validation Test

Verify that the optimization produces identical results:

```bash
Rscript scripts/validate_optimization.R
```

This runs both optimized and legacy approaches on a small dataset and compares the results.

## Troubleshooting

### Common Issues

**Issue**: `Error: Data file not found: current.csv`
- **Solution**: Ensure `current.csv` is in the project root directory, or update `config$data_file` with the correct path.

**Issue**: `Error: No observations remain after filtering`
- **Solution**: Check that `sample_start` and `sample_end` overlap with your data's date range.

**Issue**: PLS gives unexpected dimensions
- **Solution**: Enable debug mode (`config$debug <- TRUE`) to trace PLS extraction details.

**Issue**: All predictions are NA
- **Solution**: Verify that your data has enough observations before `first_forecast_idx` (default 60).

**Issue**: `Error: category_mapping_file does not exist`
- **Solution**: Create a category mapping CSV file and set `config$category_mapping_file <- "path/to/file.csv"`. See the Category Mapping section for details.

**Issue**: `Error: The following series are missing from category_mapping_file`
- **Solution**: Ensure all series in your dataset are included in the category mapping CSV file.

**Issue**: Test results contain many NA values
- **Solution**: This is expected when n_oos < 20 or when HAC variance estimates are invalid. Check that your evaluation window is large enough.

**Issue**: Workflow is taking too long
- **Solution**: Ensure you're using `compute_evaluation()` instead of separate `compute_rmse()` and `compute_tests()` calls. See the Performance Optimization section.

**Issue**: No forecasts folder in outputs
- **Solution**: Ensure `config$save_forecasts <- TRUE` (this is the default). Also verify you're using `compute_evaluation()` and passing `forecasts = evaluation$forecasts` to `save_results()`.

**Issue**: MCS results show "skipped" status
- **Solution**: This occurs when fewer than 2 methods match the M0_set specification. Check that your method patterns correctly match available methods. Use `unique(forecasts$method_id)` to see available method IDs.

**Issue**: MCS drops many rows due to NA
- **Solution**: MCS requires a balanced panel. If >5% of rows are dropped, you may have missing forecasts for some methods at certain origins. Check that all methods have forecasts at the same origins.

### Getting Help

- Check debug logs by setting `config$debug <- TRUE`
- Use `trace_origins` to inspect specific time indices
- Examine intermediate results by running components individually

## References

**Bae (2024)**: "Some variables are transformed to be stationary. Then, the transformed variables are standardized to have unit variance and mean zero. Finally, the data are screened for outliers: Any observations whose values exceed ten times the interquartile range from the median are treated as missing values. Factors are estimated only from the balanced panel with 108 predictors."

**Diebold, F. X., & Mariano, R. S. (1995)**. Comparing predictive accuracy. *Journal of Business & Economic Statistics*, 13(3), 253-263.
- Implements one-sided DM test for comparing forecast accuracy

**Clark, T. E., & West, K. D. (2007)**. Approximately normal tests for equal predictive accuracy in nested models. *Journal of Econometrics*, 138(1), 291-311.
- Implements Clark-West test accounting for nesting bias

**Newey, W. K., & West, K. D. (1987)**. A simple, positive semi-definite, heteroskedasticity and autocorrelation consistent covariance matrix. *Econometrica*, 55(3), 703-708.
- HAC standard errors with Bartlett kernel and lag selection L = h - 1

**Hansen, P. R., Lunde, A., & Nason, J. M. (2011)**. The Model Confidence Set. *Econometrica*, 79(2), 453-497.
- Model Confidence Set procedure for identifying superior forecasting models
- Stationary bootstrap for dependent data

**Bai, J., & Ng, S. (2002)**. Determining the number of factors in approximate factor models. *Econometrica*, 70(1), 191-221.
- Information criteria (IC_p1, IC_p2, IC_p3) for selecting the number of factors
- BN-BIC (IC_p3) criterion planned for data-driven factor selection

**Politis, D. N., & Romano, J. P. (1994)**. The stationary bootstrap. *Journal of the American Statistical Association*, 89(428), 1303-1313.
- Stationary bootstrap procedure used in MCS

## License

MIT License (add LICENSE file as needed)

## Contact

Sebastian Alamo Alonso - your.email@example.com
