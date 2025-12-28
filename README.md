# Factor-Based Forecasting Framework

A modular, package-like R framework for macroeconomic forecasting using factor models (PCA and PLS). Implements the methodology from Bae (2024) with support for multiple forecasting designs: DI (Direct forecast with factors), DIAR (Direct forecast with factors and AR lags), and DIAR-LAG (Direct forecast with factors and their lags).

## Features

- **Modular Architecture**: Clean separation of concerns with dedicated modules for data I/O, preprocessing, factor extraction, forecasting, and evaluation
- **Multiple Factor Methods**: Support for both PCA and PLS factor extraction
- **Multiple Forecasting Designs**: AR, DI, DIAR, DIAR-LAG models
- **Rolling and Recursive Windows**: Both estimation schemes supported
- **Configurable Workflows**: Easy-to-swap datasets and parameters via configuration system
- **Reproducible Results**: Deterministic processing with explicit seed control
- **Comprehensive Logging**: Debug and trace functionality for troubleshooting

## Project Structure

```
.
├── R/                          # Core R modules
│   ├── config.R                # Configuration system
│   ├── data_io.R               # Data loading functions
│   ├── preprocessing.R         # Data transformation and cleaning
│   ├── factors_pca.R           # PCA factor extraction
│   ├── factors_pls.R           # PLS factor extraction
│   ├── forecasting_designs.R  # Model building/fitting/prediction (AR, DI, DIAR, DIAR-LAG)
│   ├── forecasting_models.R   # High-level forecast orchestration
│   ├── evaluation.R            # RMSE computation and evaluation
│   ├── workflow.R              # Top-level workflow orchestration
│   └── utils_logging.R         # Logging utilities
├── scripts/
│   └── run_workflow.R          # Main entry point script
├── inst/extdata/               # Optional: example config files, sample data
├── tests/testthat/             # Unit and smoke tests
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
5. Compute RMSE relative to AR benchmark
6. Save results to `outputs/[timestamp]/`

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

# Run workflow
results <- run_workflow(config)

# Compute RMSE
rmse_results <- compute_rmse(results, config)

# Save results
save_results(results, rmse_results, config)

# Print summary
print_summary(results, rmse_results, config)
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

- `dataset_id`: Dataset identifier ("US_FRED" or "EURO_AREA")
- `data_file`: Path to input CSV file
- `sample_start/end`: Sample period boundaries
- `eval_start/end`: Evaluation period boundaries
- `horizons`: Forecast horizons (e.g., c(1, 6, 12, 24))
- `series_list`: Specific series to forecast (NULL = all balanced series)
- `k_max_pca/pls`: Maximum number of factors to extract
- `models`: Models to run (c("AR", "DI", "DIAR", "DIAR-LAG"))
- `schemes`: Estimation schemes (c("recursive", "rolling"))
- `factor_methods`: Factor extraction methods (c("PCA", "PLS"))
- `debug`: Enable debug logging (TRUE/FALSE)
- `trace_origins`: Time indices to trace for detailed logging

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
├── rmse_results.csv        # Full RMSE table
├── config.rds              # Configuration used for this run
├── summary.txt             # Text summary of the run
└── diarlag_h1_recursive.png  # Example plot
```

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
- `fig1_mean_rmse_k_pca.png`: Mean RMSE of k-PCA across schemes and models
- `fig2_mean_rmse_k_pls.png`: Mean RMSE of k-PLS across schemes and models

### Manual Plot Generation

Generate plots manually from existing RMSE results:

```r
# Load existing results
rmse_results <- read.csv("outputs/[run_id]/rmse_results.csv")

# Generate Bae-style figures
output_dir <- "outputs/[run_id]"
plot_bae_fig1_kpca(rmse_results, metric = "rmse_rel", save_dir = output_dir)
plot_bae_fig2_kpls(rmse_results, metric = "rmse_rel", save_dir = output_dir)
```

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

## Testing

### Smoke Test

Run a quick smoke test on a small sample:

```r
config <- config_us_default()
config$series_list <- c("INDPRO")
config$horizons <- c(1)
config$factor_methods <- c("PCA")

results <- run_workflow(config)
rmse_results <- compute_rmse(results, config)

# Verify output structure
stopifnot(nrow(rmse_results) > 0)
stopifnot(all(c("series", "h", "scheme", "model", "mse") %in% names(rmse_results)))
```

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

### Getting Help

- Check debug logs by setting `config$debug <- TRUE`
- Use `trace_origins` to inspect specific time indices
- Examine intermediate results by running components individually

## References

Bae (2024): "Some variables are transformed to be stationary. Then, the transformed variables are standardized to have unit variance and mean zero. Finally, the data are screened for outliers: Any observations whose values exceed ten times the interquartile range from the median are treated as missing values. Factors are estimated only from the balanced panel with 108 predictors."

## License

MIT License (add LICENSE file as needed)

## Contact

Sebastian Alamo Alonso - your.email@example.com
