#' Data Preprocessing Functions
#'
#' Functions for transforming, standardizing, and cleaning data.
#'
#' @name preprocessing
NULL

#' Preprocess dataset
#'
#' Applies standardization, outlier removal, and balanced panel selection.
#'
#' @param data List returned from load_dataset()
#' @param config Configuration list
#' @return A list with components:
#'   - panel_std1: data.frame with date + all standardized/cleaned predictors
#'   - panel_final: data.frame with date + balanced predictors (no NAs)
#'   - dates: vector of dates
#'   - balanced_predictors: character vector of predictor names in balanced panel
#' @export
preprocess_dataset <- function(data, config) {
  log_info("Preprocessing dataset: standardization, outlier removal, balanced panel selection", config)

  raw_df <- data$raw_data
  dates <- data$dates

  # Extract predictor matrix (all columns except date)
  X_raw <- as.matrix(raw_df[, -1])

  # Standardize: mean 0, unit variance (column-wise, NA-aware)
  if (config$standardize) {
    X_std <- scale(X_raw)
  } else {
    X_std <- X_raw
  }

  panel_std0 <- data.frame(date = dates, X_std)

  # Remove outliers: values > 10*IQR from median are set to NA
  # Bae (2024): "Any observations whose values exceed ten times the interquartile range from the median are treated as missing values"
  X_std_clean <- as.data.frame(lapply(panel_std0[, -1], function(x) {
    remove_outliers_iqr(x, multiplier = config$outlier_iqr_multiplier)
  }))

  panel_std1 <- data.frame(date = panel_std0$date, X_std_clean)

  # Create balanced panel (predictors with no NAs)
  # Bae (2024): "Factors are estimated only from the balanced panel with 108 predictors."
  if (config$require_balanced_panel) {
    predictor_matrix <- panel_std1[, -1]

    good_cols <- colnames(predictor_matrix)[
      apply(predictor_matrix, 2, function(col) all(is.finite(col)))
    ]

    log_info(sprintf("Balanced predictors: %d (out of %d total)", length(good_cols), ncol(predictor_matrix)), config)

    panel_final <- data.frame(
      date = panel_std1$date,
      predictor_matrix[, good_cols, drop = FALSE]
    )

    balanced_predictors <- good_cols
  } else {
    panel_final <- panel_std1
    balanced_predictors <- colnames(panel_std1)[-1]
  }

  list(
    panel_std1 = panel_std1,
    panel_final = panel_final,
    dates = dates,
    balanced_predictors = balanced_predictors
  )
}

#' Remove outliers based on IQR rule
#'
#' Values exceeding multiplier*IQR from the median are set to NA.
#'
#' @param x Numeric vector
#' @param multiplier IQR multiplier (default 10)
#' @return Numeric vector with outliers replaced by NA
remove_outliers_iqr <- function(x, multiplier = 10) {
  med <- median(x, na.rm = TRUE)
  iqr <- IQR(x, na.rm = TRUE)

  # If IQR is zero or NA, do nothing
  if (is.na(iqr) || iqr == 0) {
    return(x)
  }

  lower <- med - multiplier * iqr
  upper <- med + multiplier * iqr

  x[x < lower | x > upper] <- NA
  x
}

#' Construct h-step ahead targets for all series
#'
#' @param panel_final data.frame with date + balanced predictors
#' @param horizons Integer vector of forecast horizons
#' @return A nested list: targets_list[[series_name]][[paste0("h", h)]] = numeric vector of length T
#' @export
construct_targets <- function(panel_final, horizons) {
  targets_list <- list()

  for (j in 2:ncol(panel_final)) {
    series_name <- names(panel_final)[j]
    x <- panel_final[[j]]

    targets_list[[series_name]] <- lapply(horizons, function(h) {
      construct_target_h(x, h)
    })
    names(targets_list[[series_name]]) <- paste0("h", horizons)
  }

  targets_list
}

#' Construct h-step ahead target series
#'
#' For a series x of length T, the h-step target y_h[t] = x[t+h]
#'
#' @param x Numeric vector (length T)
#' @param h Forecast horizon (integer)
#' @return Numeric vector of length T, with y_h[t] = x[t+h] and NAs for t > T-h
construct_target_h <- function(x, h) {
  n <- length(x)
  out <- rep(NA_real_, n)
  out[1:(n - h)] <- x[(1 + h):n]
  out
}
