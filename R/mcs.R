#' Model Confidence Set (MCS) Evaluation
#'
#' Wrapper around the MCS package (Catania & Bernardi) which implements the
#' Hansen, Lunde, and Nason (2011) Model Confidence Set procedure for comparing
#' forecast accuracy across multiple methods.
#'
#' @name mcs
#'
#' @details
#' This module provides integration with the MCS package from CRAN.
#' The MCS package correctly implements the Hansen et al. (2011) procedure including:
#' - Pairwise loss differentials d_ij = L_i - L_j
#' - TR (range) and Tmax statistics as defined in the paper
#' - Block bootstrap with AR-based block length selection
#' - Sequential elimination procedure
#'
#' Install the MCS package with: install.packages("MCS")
#'
#' @references
#' Hansen, P. R., Lunde, A., & Nason, J. M. (2011). The Model Confidence Set.
#' Econometrica, 79(2), 453-497.
#'
#' Bernardi, M., & Catania, L. (2014). The Model Confidence Set package for R.
NULL

# ============================================================================
# Loss Functions
# ============================================================================

#' Squared Error Loss
#' @param y_true Vector of realized values
#' @param y_hat Vector of forecast values
#' @return Vector of squared errors
loss_se <- function(y_true, y_hat) {
  (y_true - y_hat)^2
}

#' Absolute Error Loss
#' @param y_true Vector of realized values
#' @param y_hat Vector of forecast values
#' @return Vector of absolute errors
loss_ae <- function(y_true, y_hat) {
  abs(y_true - y_hat)
}

#' Get loss function by name
#' @param loss_name Character: "se" or "ae"
#' @return Loss function
get_loss_fn <- function(loss_name) {

  switch(loss_name,
    "se" = loss_se,
    "ae" = loss_ae,
    stop(sprintf("Unknown loss function: %s. Use 'se' or 'ae'.", loss_name))
  )
}

# ============================================================================
# MCS Bootstrap Implementation (Hansen, Lunde, Nason 2011)
# ============================================================================

#' Stationary Bootstrap for time series
#'
#' Implements Politis & Romano (1994) stationary bootstrap.
#'
#' @param x Matrix (T x M) of data to bootstrap
#' @param B Number of bootstrap samples
#' @param expected_block_length Expected block length (geometric distribution)
#' @param seed Random seed
#' @return Array (T x M x B) of bootstrap samples
stationary_bootstrap <- function(x, B, expected_block_length = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  T_obs <- nrow(x)
  M <- ncol(x)

  # Default block length: Politis-Romano rule ~ T^(1/3)

if (is.null(expected_block_length)) {
    expected_block_length <- ceiling(T_obs^(1/3))
  }

  # Probability of starting new block
  p <- 1 / expected_block_length

  # Generate bootstrap samples
  boot_samples <- array(NA_real_, dim = c(T_obs, M, B))

  for (b in 1:B) {
    # Generate bootstrap indices
    indices <- integer(T_obs)
    indices[1] <- sample.int(T_obs, 1)

    for (t in 2:T_obs) {
      if (runif(1) < p) {
        # Start new block
        indices[t] <- sample.int(T_obs, 1)
      } else {
        # Continue current block (wrap around)
        indices[t] <- (indices[t-1] %% T_obs) + 1
      }
    }

    boot_samples[, , b] <- x[indices, , drop = FALSE]
  }

  boot_samples
}

#' Compute MCS t-statistic (TR or Tmax)
#'
#' @param d_ij Matrix (T x M) of loss differentials relative to row means
#' @param stat_type Character: "TR" (range) or "Tmax" (max)
#' @return Numeric test statistic
compute_mcs_stat <- function(d_ij, stat_type = "TR") {
  T_obs <- nrow(d_ij)
  M <- ncol(d_ij)

  # Mean loss differential for each model relative to average
  d_bar <- colMeans(d_ij)

  # Variance estimates (HAC would be better but use simple for now)
  var_d <- apply(d_ij, 2, var) / T_obs

  # t-statistics for each model
  t_stats <- d_bar / sqrt(pmax(var_d, 1e-10))

  if (stat_type == "Tmax") {
    max(t_stats)
  } else if (stat_type == "TR") {
    max(t_stats) - min(t_stats)
  } else {
    stop(sprintf("Unknown stat_type: %s", stat_type))
  }
}

#' Run MCS procedure for a single candidate set
#'
#' Wrapper around the MCS package's MCSprocedure function.
#' Falls back to a simplified internal implementation if MCS package is not installed.
#'
#' @param L Matrix (T x M) of losses, columns are methods
#' @param alpha Significance level (default: 0.10)
#' @param B Number of bootstrap replications (default: 1000)
#' @param stat_type Test statistic: "TR" or "Tmax" (default: "Tmax")
#' @param block_length Block length for bootstrap (NULL = automatic via AR fitting)
#' @param seed Random seed (not used by MCS package, included for interface compatibility)
#' @return List with:
#'   - superior_set: Character vector of methods in MCS
#'   - eliminated: Character vector of eliminated methods (in order)
#'   - pvalues: Named vector of p-values for each method
#'   - all_methods: All input methods
#'   - status: "ok", "fallback", or error status
#' @export
run_mcs <- function(L, alpha = 0.10, B = 1000, stat_type = "Tmax",
                    block_length = NULL, seed = NULL) {

  if (!is.matrix(L)) L <- as.matrix(L)

  T_obs <- nrow(L)
  M_init <- ncol(L)
  method_names <- colnames(L)

  if (is.null(method_names)) {
    method_names <- paste0("M", 1:M_init)
    colnames(L) <- method_names
  }

  if (M_init < 2) {
    return(list(
      superior_set = method_names,
      eliminated = character(0),
      pvalues = setNames(1, method_names),
      all_methods = method_names,
      status = "ok",
      message = "Only one method, trivially in MCS"
    ))
  }

  # Try to use MCS package (correct implementation)
  if (requireNamespace("MCS", quietly = TRUE)) {
    tryCatch({
      # MCS package expects Loss matrix with models as columns
      # statistic must be "Tmax" or "TR"
      mcs_result <- MCS::MCSprocedure(
        Loss = L,
        alpha = alpha,
        B = B,
        statistic = stat_type,
        k = block_length,
        verbose = FALSE
      )

      # Extract results from SSM object
      # The MCS package returns an S4 object of class "SSM" with slots:
      # @show: matrix with results (columns: Rank_M, v_M, MCS_M, Rank_R, v_R, MCS_R, Loss)
      # @Info: list with model.names, elapsed.time, statistic, n_elim, mcs_pvalue, alpha, B, k
      # @Bootstrap: list with TR and Tmax bootstrap distributions

      results_matrix <- mcs_result@show
      info <- mcs_result@Info

      # Determine which MCS column to use based on statistic
      mcs_col <- if (stat_type == "Tmax") "MCS_M" else "MCS_R"

      # Get the methods in the MCS (superior set) - MCS column = 1 means included
      if (mcs_col %in% colnames(results_matrix)) {
        in_mcs <- results_matrix[, mcs_col] == 1
        superior_set <- rownames(results_matrix)[in_mcs]
        eliminated <- rownames(results_matrix)[!in_mcs]
      } else {
        # Fallback: try to parse from first available MCS column
        mcs_cols <- grep("^MCS_", colnames(results_matrix), value = TRUE)
        if (length(mcs_cols) > 0) {
          in_mcs <- results_matrix[, mcs_cols[1]] == 1
          superior_set <- rownames(results_matrix)[in_mcs]
          eliminated <- rownames(results_matrix)[!in_mcs]
        } else {
          stop("Cannot find MCS column in results matrix")
        }
      }

      # Extract p-value from Info slot
      mcs_pvalue <- if (!is.null(info$mcs_pvalue)) info$mcs_pvalue else NA_real_

      # Create p-values vector (MCS package gives one overall p-value, not per-method)
      # Methods in MCS get p-value = mcs_pvalue, eliminated get sequential p-values
      # For simplicity, assign mcs_pvalue to superior set, NA to eliminated
      pvalues <- setNames(rep(NA_real_, length(method_names)), method_names)
      pvalues[superior_set] <- mcs_pvalue

      return(list(
        superior_set = superior_set,
        eliminated = eliminated,
        pvalues = pvalues,
        all_methods = method_names,
        alpha = alpha,
        stat_type = stat_type,
        B = B,
        status = "ok",
        message = sprintf("MCS contains %d of %d methods (via MCS package, p=%.4f)",
                          length(superior_set), M_init, mcs_pvalue)
      ))
    }, error = function(e) {
      warning(sprintf("MCS package failed: %s. Using fallback.", e$message))
    })
  } else {
    warning("MCS package not installed. Using simplified fallback implementation. ",
            "For correct results, install with: install.packages('MCS')")
  }

  # Fallback: simplified implementation (may not match paper exactly)
  # This is kept for environments where MCS package cannot be installed
  run_mcs_fallback(L, alpha, B, stat_type, block_length, seed)
}

#' Fallback MCS implementation (simplified)
#'
#' A simplified MCS implementation used when the MCS package is not available.
#' WARNING: This may not exactly match the Hansen et al. (2011) procedure.
#' For correct results, install the MCS package.
#'
#' @param L Matrix (T x M) of losses
#' @param alpha Significance level
#' @param B Bootstrap replications
#' @param stat_type "TR" or "Tmax"
#' @param block_length Block length (NULL = T^(1/3))
#' @param seed Random seed
#' @return List with MCS results
#' @keywords internal
run_mcs_fallback <- function(L, alpha = 0.10, B = 1000, stat_type = "Tmax",
                              block_length = NULL, seed = NULL) {

  T_obs <- nrow(L)
  M_init <- ncol(L)
  method_names <- colnames(L)

  # Initialize
 current_set <- method_names
  eliminated <- character(0)
  pvalues <- setNames(rep(NA_real_, M_init), method_names)

  while (length(current_set) > 1) {
    # Subset loss matrix to current set
    L_current <- L[, current_set, drop = FALSE]
    M_current <- length(current_set)

    # Compute loss differentials relative to row mean
    row_means <- rowMeans(L_current)
    d_ij <- L_current - row_means

    # Observed test statistic
    stat_obs <- compute_mcs_stat(d_ij, stat_type)

    # Bootstrap distribution
    if (!is.null(seed)) set.seed(seed + length(eliminated))
    boot_samples <- stationary_bootstrap(d_ij, B, block_length, seed = NULL)

    boot_stats <- numeric(B)
    for (b in 1:B) {
      boot_stats[b] <- compute_mcs_stat(boot_samples[, , b], stat_type)
    }

    # P-value (proportion of bootstrap stats >= observed)
    pval <- mean(boot_stats >= stat_obs)

    if (pval > alpha) {
      # Cannot reject null: current set is the MCS
      break
    }

    # Find and eliminate worst model (highest mean loss differential)
    mean_d <- colMeans(d_ij)
    worst_idx <- which.max(mean_d)
    worst_method <- current_set[worst_idx]

    pvalues[worst_method] <- pval
    eliminated <- c(eliminated, worst_method)
    current_set <- setdiff(current_set, worst_method)
  }

  # Assign p-value = 1 to methods in superior set (not eliminated)
  pvalues[current_set] <- 1

  list(
    superior_set = current_set,
    eliminated = eliminated,
    pvalues = pvalues,
    all_methods = method_names,
    alpha = alpha,
    stat_type = stat_type,
    B = B,
    status = "fallback",
    message = sprintf("MCS contains %d of %d methods (FALLBACK - install MCS package for correct results)",
                      length(current_set), M_init)
  )
}

# ============================================================================
# Method ID Resolution
# ============================================================================

#' Resolve method specification to actual method_ids
#'
#' Converts symbolic method names (like "k1-PLS", "AR", "PCA_BNBIC", "PLS_ON")
#' to actual method_id patterns that can be matched against the forecasts data.
#'
#' Supported formats:
#' - New grid: {scheme}_{factor_spec_id}_{model_class}_k{n} (e.g., recursive_PCA_grid_DI_k3)
#' - New dynamic: {scheme}_{factor_spec_id}_{model_class} (e.g., recursive_PCA_BNBIC_DI)
#'
#' @param method_spec Character vector of method specifications
#' @param available_methods Character vector of available method_ids
#' @param scheme Scheme filter ("recursive" or "rolling")
#' @param model_class Model class filter (e.g., "DI", "DIAR")
#' @return Character vector of matched method_ids
#' @export
resolve_method_specs <- function(method_spec, available_methods, scheme = NULL, model_class = NULL) {
  resolved <- character(0)

  for (spec in method_spec) {
    pattern <- NULL

    # Handle symbolic names
    if (spec == "AR") {
      if (!is.null(scheme)) {
        pattern <- paste0("^", scheme, "_AR$")
      } else {
        pattern <- "_AR$"
      }

    } else if (grepl("^k([0-9]+)-PLS$", spec)) {
      k <- sub("^k([0-9]+)-PLS$", "\\1", spec)
      base <- if (!is.null(scheme)) paste0("^", scheme, "_PLS") else "_PLS"
      if (!is.null(model_class)) {
        pattern <- paste0(base, "_grid_", model_class, "_k", k, "$")
      } else {
        pattern <- paste0(base, "_grid_.*_k", k, "$")
      }

    } else if (grepl("^k([0-9]+)-PCA$", spec)) {
      k <- sub("^k([0-9]+)-PCA$", "\\1", spec)
      base <- if (!is.null(scheme)) paste0("^", scheme, "_PCA") else "_PCA"
      if (!is.null(model_class)) {
        pattern <- paste0(base, "_grid_", model_class, "_k", k, "$")
      } else {
        pattern <- paste0(base, "_grid_.*_k", k, "$")
      }

    } else if (spec == "PCA_BNBIC" || spec == "PCA-BNBIC") {
      # Dynamic PCA with BN-BIC
      base <- if (!is.null(scheme)) paste0("^", scheme, "_PCA_BNBIC_") else "_PCA_BNBIC_"
      if (!is.null(model_class)) {
        pattern <- paste0(base, model_class, "$")
      } else {
        pattern <- paste0(base, ".*$")
      }

    } else if (spec == "PLS_ON" || spec == "PLS-ON") {
      # Dynamic PLS with Onatski
      base <- if (!is.null(scheme)) paste0("^", scheme, "_PLS_ON_") else "_PLS_ON_"
      if (!is.null(model_class)) {
        pattern <- paste0(base, model_class, "$")
      } else {
        pattern <- paste0(base, ".*$")
      }

    } else if (spec == "PCA_grid") {
      # All grid PCA methods
      base <- if (!is.null(scheme)) paste0("^", scheme, "_PCA_grid_") else "_PCA_grid_"
      if (!is.null(model_class)) {
        pattern <- paste0(base, model_class, "_k[0-9]+$")
      } else {
        pattern <- paste0(base, ".*_k[0-9]+$")
      }

    } else if (spec == "PLS_grid") {
      # All grid PLS methods
      base <- if (!is.null(scheme)) paste0("^", scheme, "_PLS_grid_") else "_PLS_grid_"
      if (!is.null(model_class)) {
        pattern <- paste0(base, model_class, "_k[0-9]+$")
      } else {
        pattern <- paste0(base, ".*_k[0-9]+$")
      }

    } else {
      # Treat as regex pattern or literal match
      pattern <- spec
    }

    if (!is.null(pattern)) {
      matches <- grep(pattern, available_methods, value = TRUE)
      resolved <- c(resolved, matches)
    }
  }

  unique(resolved)
}

# ============================================================================
# Main MCS Evaluation Function
# ============================================================================

#' Run MCS evaluation across all configured slices
#'
#' Runs MCS for each slice (series_id, h, scheme, model_class) and each M0_set.
#' Slices by model_class (forecast equation) so that DI models are compared with DI,
#' DIAR with DIAR, etc. If config$mcs$include_ar_in_M0 = TRUE (default), the AR
#' benchmark is automatically included in each MCS comparison.
#'
#' @param forecasts Tibble of forecasts (from compute_evaluation)
#' @param config Configuration list with mcs settings
#' @return Tibble with MCS results for all slices and M0_sets
#' @export
#' @importFrom dplyr filter group_by summarise mutate select distinct bind_rows
#' @importFrom tidyr pivot_wider
compute_mcs_evaluation <- function(forecasts, config) {

  if (!isTRUE(config$mcs$enabled)) {
    return(NULL)
  }

  log_info("Computing MCS evaluation", config)

  mcs_cfg <- config$mcs
  loss_fn <- get_loss_fn(mcs_cfg$loss)
  include_ar <- isTRUE(mcs_cfg$include_ar_in_M0)

  # Get unique slices - slice by model_class (forecast equation)
  # We exclude AR from slicing since it will be added separately if include_ar_in_M0 = TRUE
  factor_forecasts <- forecasts %>%
    dplyr::filter(model_class != "AR")

  slices <- factor_forecasts %>%
    dplyr::select(run_id, series_id, h, scheme, model_class) %>%
    dplyr::distinct()

  log_info(sprintf("MCS evaluation: %d slices, %d M0_sets, %d alphas, include_ar=%s",
                   nrow(slices), length(mcs_cfg$M0_sets), length(mcs_cfg$alphas),
                   if (include_ar) "TRUE" else "FALSE"), config)

  results_list <- list()

  for (i in 1:nrow(slices)) {
    slice <- slices[i, ]

    # Filter forecasts for this slice (specific model_class)
    slice_data <- forecasts %>%
      dplyr::filter(
        run_id == slice$run_id,
        series_id == slice$series_id,
        h == slice$h,
        scheme == slice$scheme,
        model_class == slice$model_class
      )

    # If include_ar_in_M0, also get AR forecasts for this slice
    if (include_ar) {
      ar_data <- forecasts %>%
        dplyr::filter(
          run_id == slice$run_id,
          series_id == slice$series_id,
          h == slice$h,
          scheme == slice$scheme,
          model_class == "AR"
        )
      slice_data <- dplyr::bind_rows(slice_data, ar_data)
    }

    available_methods <- unique(slice_data$method_id)

    if (length(available_methods) < 2) {
      # Skip slice with fewer than 2 methods
      for (m0_name in names(mcs_cfg$M0_sets)) {
        for (alpha in mcs_cfg$alphas) {
          results_list[[length(results_list) + 1]] <- tibble::tibble(
            run_id = slice$run_id,
            series_id = slice$series_id,
            h = slice$h,
            scheme = slice$scheme,
            model_class = slice$model_class,
            alpha = alpha,
            M0_set_id = m0_name,
            loss = mcs_cfg$loss,
            test_stat = mcs_cfg$test_stat,
            B = mcs_cfg$B,
            n_methods_input = length(available_methods),
            n_methods_mcs = NA_integer_,
            included_methods = paste(available_methods, collapse = ";"),
            superior_set = NA_character_,
            pvalues = NA_character_,
            status = "skipped",
            message = "Fewer than 2 methods available"
          )
        }
      }
      next
    }

    # Process each M0_set
    for (m0_name in names(mcs_cfg$M0_sets)) {
      m0_spec <- mcs_cfg$M0_sets[[m0_name]]

      # Resolve method specs to actual method_ids
      # Filter by model_class for factor-based methods
      resolved_methods <- resolve_method_specs(
        m0_spec,
        available_methods,
        scheme = slice$scheme,
        model_class = slice$model_class
      )

      # If include_ar_in_M0 = TRUE, ensure AR is included
      if (include_ar) {
        ar_method_id <- paste0(slice$scheme, "_AR")
        if (ar_method_id %in% available_methods && !(ar_method_id %in% resolved_methods)) {
          resolved_methods <- c(resolved_methods, ar_method_id)
        }
      }

      if (length(resolved_methods) < 2) {
        for (alpha in mcs_cfg$alphas) {
          results_list[[length(results_list) + 1]] <- tibble::tibble(
            run_id = slice$run_id,
            series_id = slice$series_id,
            h = slice$h,
            scheme = slice$scheme,
            model_class = slice$model_class,
            alpha = alpha,
            M0_set_id = m0_name,
            loss = mcs_cfg$loss,
            test_stat = mcs_cfg$test_stat,
            B = mcs_cfg$B,
            n_methods_input = length(resolved_methods),
            n_methods_mcs = NA_integer_,
            included_methods = paste(resolved_methods, collapse = ";"),
            superior_set = NA_character_,
            pvalues = NA_character_,
            status = "skipped",
            message = sprintf("M0_set '%s' resolved to fewer than 2 methods", m0_name)
          )
        }
        next
      }

      # Build loss matrix for this slice and M0_set
      tryCatch({
        # Filter to resolved methods and compute losses
        L_data <- slice_data %>%
          dplyr::filter(method_id %in% resolved_methods) %>%
          dplyr::filter(!is.na(y_true) & !is.na(y_hat)) %>%
          dplyr::mutate(loss_val = loss_fn(y_true, y_hat)) %>%
          dplyr::select(origin_index, method_id, loss_val)

        # Pivot to matrix
        L_wide <- L_data %>%
          tidyr::pivot_wider(
            id_cols = origin_index,
            names_from = method_id,
            values_from = loss_val
          ) %>%
          dplyr::arrange(origin_index)

        # Complete cases only
        L_matrix <- as.matrix(L_wide[, -1, drop = FALSE])
        complete_mask <- stats::complete.cases(L_matrix)
        L_matrix <- L_matrix[complete_mask, , drop = FALSE]

        if (nrow(L_matrix) < 20) {
          for (alpha in mcs_cfg$alphas) {
            results_list[[length(results_list) + 1]] <- tibble::tibble(
              run_id = slice$run_id,
              series_id = slice$series_id,
              h = slice$h,
              scheme = slice$scheme,
              model_class = slice$model_class,
              alpha = alpha,
              M0_set_id = m0_name,
              loss = mcs_cfg$loss,
              test_stat = mcs_cfg$test_stat,
              B = mcs_cfg$B,
              n_methods_input = ncol(L_matrix),
              n_methods_mcs = NA_integer_,
              included_methods = paste(colnames(L_matrix), collapse = ";"),
              superior_set = NA_character_,
              pvalues = NA_character_,
              status = "skipped",
              message = sprintf("Insufficient observations after alignment: %d", nrow(L_matrix))
            )
          }
          next
        }

        # Run MCS for each alpha
        for (alpha in mcs_cfg$alphas) {
          mcs_result <- run_mcs(
            L = L_matrix,
            alpha = alpha,
            B = mcs_cfg$B,
            stat_type = mcs_cfg$test_stat,
            block_length = mcs_cfg$block_length,
            seed = mcs_cfg$seed
          )

          results_list[[length(results_list) + 1]] <- tibble::tibble(
            run_id = slice$run_id,
            series_id = slice$series_id,
            h = slice$h,
            scheme = slice$scheme,
            model_class = slice$model_class,
            alpha = alpha,
            M0_set_id = m0_name,
            loss = mcs_cfg$loss,
            test_stat = mcs_cfg$test_stat,
            B = mcs_cfg$B,
            n_methods_input = length(mcs_result$all_methods),
            n_methods_mcs = length(mcs_result$superior_set),
            included_methods = paste(mcs_result$all_methods, collapse = ";"),
            superior_set = paste(mcs_result$superior_set, collapse = ";"),
            pvalues = paste(sprintf("%s=%.4f", names(mcs_result$pvalues), mcs_result$pvalues), collapse = ";"),
            status = mcs_result$status,
            message = mcs_result$message
          )
        }

      }, error = function(e) {
        for (alpha in mcs_cfg$alphas) {
          results_list[[length(results_list) + 1]] <<- tibble::tibble(
            run_id = slice$run_id,
            series_id = slice$series_id,
            h = slice$h,
            scheme = slice$scheme,
            model_class = slice$model_class,
            alpha = alpha,
            M0_set_id = m0_name,
            loss = mcs_cfg$loss,
            test_stat = mcs_cfg$test_stat,
            B = mcs_cfg$B,
            n_methods_input = NA_integer_,
            n_methods_mcs = NA_integer_,
            included_methods = NA_character_,
            superior_set = NA_character_,
            pvalues = NA_character_,
            status = "error",
            message = conditionMessage(e)
          )
        }
      })
    }
  }

  mcs_results <- dplyr::bind_rows(results_list)

  n_ok <- sum(mcs_results$status == "ok")
  n_fallback <- sum(mcs_results$status == "fallback")
  n_skipped <- sum(mcs_results$status == "skipped")
  n_error <- sum(mcs_results$status == "error")

  log_info(sprintf("MCS evaluation complete: %d results (%d ok, %d fallback, %d skipped, %d errors)",
                   nrow(mcs_results), n_ok, n_fallback, n_skipped, n_error), config)

  if (n_fallback > 0) {
    log_warn(sprintf("MCS fallback used %d times. Install MCS package for correct results: install.packages('MCS')",
                     n_fallback), config)
  }

  mcs_results
}
