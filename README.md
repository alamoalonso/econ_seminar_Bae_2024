# Factor-Based Forecasting Framework

A modular R framework for macroeconomic forecasting using factor models (PCA and PLS). Implements the methodology from Bae (2024) with support for multiple forecasting designs: DI, DIAR, and DIAR-LAG.

## Features

- **Multiple Factor Methods**: PCA and PLS factor extraction
- **Multiple Forecasting Designs**: AR, DI, DIAR, DIAR-LAG models
- **Rolling and Recursive Windows**: Both estimation schemes supported
- **Forecast Comparison Tests**: Diebold-Mariano and Clark-West tests with HAC standard errors
- **Model Confidence Set (MCS)**: Hansen-Lunde-Nason (2011) procedure
- **Configurable Workflows**: Easy-to-swap datasets and parameters via configuration system
- **Reproducible Results**: Deterministic processing with explicit seed control

## Project Structure

```
.
├── R/                          # Core R modules
│   ├── config.R                # Configuration system
│   ├── data_io.R               # Data loading functions
│   ├── preprocessing.R         # Data transformation and cleaning
│   ├── factors_pca.R           # PCA factor extraction
│   ├── factors_pls.R           # PLS factor extraction
│   ├── paper_factors_pca.R     # Paper-compliant PCA (Bae 2024 normalization)
│   ├── paper_factors_pls.R     # Paper-compliant PLS
│   ├── factor_specs.R          # Factor spec creation and validation
│   ├── k_selection.R           # BN-BIC and Onatski k-selection rules
│   ├── factor_cache.R          # Per-origin caching for efficient extraction
│   ├── forecasting_designs.R   # Model building/fitting/prediction (AR, DI, DIAR, DIAR-LAG)
│   ├── forecasting_models.R    # High-level forecast orchestration
│   ├── evaluation.R            # RMSE, DM/CW tests, unified evaluation
│   ├── mcs.R                   # Model Confidence Set implementation
│   ├── plots_rmse.R            # Bae-style RMSE plotting functions
│   ├── workflow.R              # Top-level workflow orchestration
│   └── utils_logging.R         # Logging utilities
├── scripts/
│   └── run_workflow.R          # Main entry point script
├── tests/testthat/             # Unit and integration tests
├── outputs/                    # Results output directory (created automatically)
├── DESCRIPTION                 # Package metadata
├── NAMESPACE                   # Package namespace
└── README.md                   # This file
```

## Installation

This project is structured as an R package but does not need to be formally installed. Ensure all dependencies are available:

```r
install.packages(c(
  "fbi",        # For FRED-MD data handling
  "pls",        # For PLS factor extraction
  "dplyr",
  "tibble",
  "readr",
  "ggplot2",
  "patchwork",  # For combining plots
  "purrr",
  "lubridate",
  "readxl",
  "janitor"
))

# Recommended: For Model Confidence Set evaluation
install.packages("MCS")
```

## Quick Start

### Running the US Workflow

```bash
# From the project root directory
Rscript scripts/run_workflow.R US
```

### Running from R Console

```r
# Set working directory to project root
setwd("path/to/project")

# Source all modules
source_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (f in source_files) source(f)

library(fbi); library(pls); library(tidyverse); library(lubridate)

# Create configuration and run
config <- config_us_default()
results <- run_workflow(config)

# Compute evaluation (RMSE + DM/CW tests)
evaluation <- compute_evaluation(results, config)
rmse_results  <- evaluation$rmse_results
tests_results <- evaluation$tests_results
forecasts     <- evaluation$forecasts

# Optional: MCS
if (isTRUE(config$mcs$enabled) && !is.null(forecasts)) {
  mcs_results <- compute_mcs_evaluation(forecasts, config)
}

# Save and summarize
save_results(results, rmse_results, config, tests_results,
             forecasts = forecasts, mcs_results = mcs_results)
print_summary(results, rmse_results, config, tests_results)
```

## Configuration

Two default configurations are provided:

```r
config <- config_us_default()   # US FRED-MD (1959-03 to 2019-12)
config <- config_euro_default() # Euro Area (stub — update data_file and dates)
```

Key parameters:
- `horizons`: Forecast horizons, e.g. `c(1, 6, 12, 24)`
- `k_max_pca` / `k_max_pls`: Maximum factors to extract
- `models`: Models to run, e.g. `c("AR", "DI", "DIAR", "DIAR-LAG")`
- `schemes`: Estimation schemes `c("recursive", "rolling")`
- `do_tests`: Enable DM/CW tests (default: `TRUE`)
- `mcs$enabled`: Enable MCS evaluation (default: `FALSE`)
- `make_plots`: Generate Bae-style figures (default: `TRUE`)
- `save_forecasts`: Persist OOS forecasts for MCS (default: `TRUE`)

## Output Structure

Results are saved to `outputs/[run_id]/`:

```
outputs/20231214_153045/
├── rmse_results.csv          # Full RMSE table
├── tests_results.csv         # DM/CW test results
├── table5_dm.csv             # Table 5 summary (Diebold-Mariano)
├── table5_cw.csv             # Table 5 summary (Clark-West)
├── forecasts/                # OOS forecasts for MCS (partitioned parquet)
├── config.rds                # Configuration used for this run
├── summary.txt               # Text summary
├── fig1_mean_rmse_k_pca.png  # PCA RMSE plot
└── fig2_mean_rmse_k_pls.png  # PLS RMSE plot
```

## References

**Bae, J. (2024).** *Forecasting with factor models using PLS*. Working paper.

**Hansen, P. R., Lunde, A., & Nason, J. M. (2011).** The Model Confidence Set. *Econometrica*, 79(2), 453–497.

**Bai, J., & Ng, S. (2002).** Determining the number of factors in approximate factor models. *Econometrica*, 70(1), 191–221.
