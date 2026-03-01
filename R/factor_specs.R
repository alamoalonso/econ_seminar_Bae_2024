#' Factor Specification Utilities
#'
#' Functions for creating, validating, and managing factor extraction specifications.
#' Factor specs enable mixed grid + dynamic k-selection in a single workflow run.
#'
#' @name factor_specs
NULL

#' Create a factor specification
#'
#' @param id Unique identifier for this spec (used in method_id)
#' @param factor_method Factor extraction method: "PCA", "PLS", or "1-PLS"
#' @param k_mode Mode of k selection: "grid" (evaluate k=1..k_max) or "dynamic" (select k_hat per origin)
#' @param k_rule Decision rule for dynamic mode: NULL for grid, "bn_bic" or "onatski" for dynamic
#' @param k_max Maximum number of factors (upper bound for both grid and dynamic)
#' @return A validated factor_spec list
#' @export
#' @examples
#' # Grid spec: evaluate k = 1..12
#' spec_grid <- create_factor_spec("PCA_grid", "PCA", "grid", NULL, 12)
#'
#' # Dynamic spec: BN-BIC selection
#' spec_bnbic <- create_factor_spec("PCA_BNBIC", "PCA", "dynamic", "bn_bic", 12)
create_factor_spec <- function(id, factor_method, k_mode, k_rule = NULL, k_max) {
  # Validate factor_method

valid_factor_methods <- c("PCA", "PLS", "1-PLS")
  if (!factor_method %in% valid_factor_methods) {
    stop(sprintf("Invalid factor_method '%s'. Must be one of: %s",
                 factor_method, paste(valid_factor_methods, collapse = ", ")))
  }

  # Validate k_mode
  valid_k_modes <- c("grid", "dynamic")
  if (!k_mode %in% valid_k_modes) {
    stop(sprintf("Invalid k_mode '%s'. Must be one of: %s",
                 k_mode, paste(valid_k_modes, collapse = ", ")))
  }

  # Validate k_rule based on k_mode
  if (k_mode == "grid") {
    if (!is.null(k_rule)) {
      stop("k_rule must be NULL for k_mode = 'grid'")
    }
  } else {
    # k_mode == "dynamic"
    valid_k_rules <- c("bn_bic", "onatski")
    if (is.null(k_rule) || !k_rule %in% valid_k_rules) {
      stop(sprintf("k_rule must be one of '%s' for k_mode = 'dynamic'",
                   paste(valid_k_rules, collapse = "', '")))
    }

    # Validate rule/method compatibility
    if (k_rule == "bn_bic" && factor_method != "PCA") {
      stop("bn_bic rule is only valid for PCA factor_method")
    }
    if (k_rule == "onatski" && !factor_method %in% c("PLS", "1-PLS")) {
      stop("onatski rule is only valid for PLS or 1-PLS factor_method")
    }
  }

  # Validate k_max
  if (!is.numeric(k_max) || length(k_max) != 1 || k_max < 1 || k_max != as.integer(k_max)) {
    stop("k_max must be a positive integer >= 1")
  }

  # Validate id
  if (!is.character(id) || length(id) != 1 || nchar(id) == 0) {
    stop("id must be a non-empty character string")
  }

  list(
    id = id,
    factor_method = factor_method,
    k_mode = k_mode,
    k_rule = k_rule,
    k_max = as.integer(k_max)
  )
}

#' Expand factor_methods to factor_specs (backward compatibility)
#'
#' Converts the simple factor_methods vector to a list of grid factor_specs.
#' This ensures backward compatibility for configs without explicit factor_specs.
#'
#' @param factor_methods Character vector of factor methods (e.g., c("PCA", "PLS"))
#' @param k_max_pca Integer: maximum k for PCA
#' @param k_max_pls Integer: maximum k for PLS
#' @return List of factor_spec objects (grid mode)
#' @export
expand_factor_methods_to_specs <- function(factor_methods, k_max_pca, k_max_pls) {
  specs <- lapply(factor_methods, function(fm) {
    k_max <- if (fm == "PCA") k_max_pca else k_max_pls
    create_factor_spec(
      id = paste0(fm, "_grid"),
      factor_method = fm,
      k_mode = "grid",
      k_rule = NULL,
      k_max = k_max
    )
  })

  names(specs) <- sapply(specs, function(s) s$id)
  specs
}

#' Validate a list of factor_specs
#'
#' Checks that all specs are valid and have unique ids.
#'
#' @param factor_specs List of factor_spec objects
#' @return TRUE if valid, stops with error otherwise
#' @export
validate_factor_specs <- function(factor_specs) {
  if (!is.list(factor_specs) || length(factor_specs) == 0) {
    stop("factor_specs must be a non-empty list")
  }

  # Check each spec has required fields
  required_fields <- c("id", "factor_method", "k_mode", "k_max")
  for (i in seq_along(factor_specs)) {
    spec <- factor_specs[[i]]
    if (!is.list(spec)) {
      stop(sprintf("factor_specs[[%d]] is not a list", i))
    }
    missing <- setdiff(required_fields, names(spec))
    if (length(missing) > 0) {
      stop(sprintf("factor_specs[[%d]] is missing required fields: %s",
                   i, paste(missing, collapse = ", ")))
    }
  }

  # Check for duplicate ids
  ids <- sapply(factor_specs, function(s) s$id)
  if (anyDuplicated(ids)) {
    dup_ids <- ids[duplicated(ids)]
    stop(sprintf("Duplicate factor_spec ids: %s", paste(unique(dup_ids), collapse = ", ")))
  }

  invisible(TRUE)
}

#' Get factor_specs from config (with backward compatibility)
#'
#' Returns factor_specs from config, auto-generating from factor_methods if needed.
#'
#' @param config Configuration list
#' @return List of factor_spec objects
#' @export
get_factor_specs <- function(config) {
  if (!is.null(config$factor_specs)) {
    # Explicit factor_specs provided
    validate_factor_specs(config$factor_specs)
    return(config$factor_specs)
  }

  # Backward compatibility: expand factor_methods to grid specs
  if (is.null(config$factor_methods)) {
    stop("Neither factor_specs nor factor_methods is defined in config")
  }

  expand_factor_methods_to_specs(
    config$factor_methods,
    config$k_max_pca,
    config$k_max_pls
  )
}

#' Get maximum k_max for a factor_method across all specs
#'
#' Used to determine how many factors to extract at each origin.
#'
#' @param factor_specs List of factor_spec objects
#' @param factor_method Character: "PCA", "PLS", or "1-PLS"
#' @return Integer: maximum k_max across matching specs, or 0 if none
#' @export
get_max_k_for_method <- function(factor_specs, factor_method) {
  matching <- Filter(function(s) s$factor_method == factor_method, factor_specs)
  if (length(matching) == 0) return(0L)
  max(sapply(matching, function(s) s$k_max))
}

#' Check if any spec uses a dynamic rule
#'
#' @param factor_specs List of factor_spec objects
#' @param rule Character: "bn_bic" or "onatski"
#' @return Logical
#' @export
has_dynamic_rule <- function(factor_specs, rule) {
  any(sapply(factor_specs, function(s) {
    s$k_mode == "dynamic" && identical(s$k_rule, rule)
  }))
}

#' Get specs for a specific factor_method
#'
#' @param factor_specs List of factor_spec objects
#' @param factor_method Character: "PCA", "PLS", or "1-PLS"
#' @return Filtered list of specs
#' @export
filter_specs_by_method <- function(factor_specs, factor_method) {
  Filter(function(s) s$factor_method == factor_method, factor_specs)
}

#' Print factor_spec summary
#'
#' @param factor_specs List of factor_spec objects
#' @export
print_factor_specs <- function(factor_specs) {
  cat("Factor Specifications:\n")
  cat(sprintf("  Total specs: %d\n", length(factor_specs)))

  for (spec in factor_specs) {
    if (spec$k_mode == "grid") {
      cat(sprintf("  - %s: %s, grid k=1..%d\n",
                  spec$id, spec$factor_method, spec$k_max))
    } else {
      cat(sprintf("  - %s: %s, dynamic (%s), k_max=%d\n",
                  spec$id, spec$factor_method, spec$k_rule, spec$k_max))
    }
  }
}
