#!/usr/bin/env Rscript

# ============================================================================
# Table: Frequency of k Attaining Minimum OOS RMSE  (PLS only)
# ============================================================================
# Generates one table per horizon x scheme combination.
#
# Input:  evaluation$rmse_results  (must exist in the R environment)
# Output: outputs/tables/kmin_counts_pls_h<h>_<scheme>.{csv,tex}
#
# Table layout (one table per h x scheme):
#   Rows   : DI, DIAR, DIAR-LAG
#   Columns: k = 1, ..., 12
#   Values : count of targets for which k minimises relative OOS RMSE
#            DIAR-LAG k > 4 shown as "--" (inadmissible)
#
# Caption in LaTeX: Factor method: PLS, Horizon: h = <h>, Scheme: <scheme>
# ============================================================================

if (basename(getwd()) == "scripts") setwd("..")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(knitr)
  library(kableExtra)
  library(readr)
})

# ============================================================================
# Parameters
# ============================================================================

HORIZONS      <- c(1, 6, 12)
K_MAX         <- 12
K_MAX_DIARLAG <- 4
MODELS        <- c("DI", "DIAR", "DIAR-LAG")
SCHEMES       <- c("recursive", "rolling")

out_dir <- file.path("outputs", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# Input validation
# ============================================================================

if (!exists("evaluation") || is.null(evaluation[["rmse_results"]])) {
  stop(
    "evaluation$rmse_results not found.\n",
    "Run the workflow first: evaluation <- compute_evaluation(results, config)"
  )
}

rmse_results <- evaluation$rmse_results

required_cols <- c("series", "h", "scheme", "factor_method", "model", "k", "rmse_rel")
missing_cols  <- setdiff(required_cols, names(rmse_results))
if (length(missing_cols) > 0) {
  stop("rmse_results is missing columns: ", paste(missing_cols, collapse = ", "))
}

available_horizons <- intersect(HORIZONS, unique(rmse_results$h))
if (length(available_horizons) == 0) {
  stop("None of the requested horizons found. Available: ",
       paste(sort(unique(rmse_results$h)), collapse = ", "))
}
HORIZONS <- sort(available_horizons)

cat("Targets :", length(unique(rmse_results$series)), "\n")
cat("Horizons:", paste(HORIZONS, collapse = ", "), "\n\n")

# ============================================================================
# Core function: build one wide table for a single (h, scheme)
# ============================================================================

build_table <- function(h_val, sch) {

  d <- rmse_results |>
    filter(
      factor_method == "PLS",
      scheme        == sch,
      model         %in% MODELS,
      h             == h_val,
      !is.na(k),
      is.finite(rmse_rel)
    )

  if (nrow(d) == 0) return(NULL)

  # For each (series, model): find the k that minimises rmse_rel.
  # Ties broken by smallest k.
  kmin <- d |>
    group_by(series, model) |>
    slice_min(rmse_rel, n = 1, with_ties = FALSE) |>
    ungroup() |>
    select(series, model, k_opt = k)

  # Count how many targets chose each k, per model
  counts <- kmin |>
    group_by(model, k_opt) |>
    summarise(count = n(), .groups = "drop")

  n_targets <- n_distinct(kmin$series)

  # Complete grid of valid (model, k) pairs — fills zeros for unoccupied k
  valid_grid <- bind_rows(
    expand_grid(model = c("DI", "DIAR"), k_opt = seq_len(K_MAX)),
    expand_grid(model = "DIAR-LAG",      k_opt = seq_len(K_MAX_DIARLAG))
  )

  # Wide table: rows = model, columns = k
  wide <- valid_grid |>
    left_join(counts, by = c("model", "k_opt")) |>
    mutate(count = replace_na(count, 0L)) |>
    pivot_wider(names_from = k_opt, values_from = count) |>
    # k > K_MAX_DIARLAG for DIAR-LAG are absent from the grid -> NA after pivot
    mutate(model = factor(model, levels = MODELS)) |>
    arrange(model)

  # Ensure all k=1..K_MAX columns exist (add NA for any missing)
  missing_k <- setdiff(as.character(seq_len(K_MAX)), names(wide))
  for (col in missing_k) wide[[col]] <- NA_integer_

  wide <- wide |>
    select(model, all_of(as.character(seq_len(K_MAX))))

  list(wide = wide, n_targets = n_targets)
}

# ============================================================================
# LaTeX formatter for one table
# ============================================================================

make_latex <- function(wide, h_val, sch, n_targets) {
  k_cols       <- as.character(seq_len(K_MAX))
  scheme_label <- if (sch == "recursive") "Recursive" else "Rolling"

  # Replace NA (inadmissible k) with "--"; keep 0 as character "0"
  display <- wide |>
    mutate(
      model = as.character(model),
      across(all_of(k_cols), ~ ifelse(is.na(.), "--", as.character(.)))
    )

  col_names <- c("Equation", paste0("$", seq_len(K_MAX), "$"))

  caption <- sprintf(
    paste0(
      "Number of targets (out of $N = %d$) for which $k$ minimises ",
      "the out-of-sample relative RMSE. ",
      "Factor method: PLS. Horizon: $h = %d$. Scheme: %s. ",
      "DIAR-LAG is estimated with at most $k = %d$ factors; ",
      "dashes (\\text{--}) indicate inadmissible values of $k$."
    ),
    n_targets, h_val, scheme_label, K_MAX_DIARLAG
  )

  kbl(
    display,
    format    = "latex",
    caption   = caption,
    col.names = col_names,
    booktabs  = TRUE,
    escape    = FALSE,
    align     = c("l", rep("r", K_MAX))
  ) |>
    kable_styling(latex_options = c("hold_position", "scale_down")) |>
    add_header_above(
      c(" " = 1, setNames(K_MAX, "Number of factors $k$")),
      escape = FALSE
    )
}

# ============================================================================
# Main loop: one table per h x scheme
# ============================================================================

for (h_val in HORIZONS) {
  for (sch in SCHEMES) {

    result <- build_table(h_val, sch)

    if (is.null(result)) {
      cat(sprintf("h=%d, %s: no data, skipped\n", h_val, sch))
      next
    }

    wide      <- result$wide
    n_targets <- result$n_targets
    stem      <- sprintf("kmin_counts_pls_h%d_%s", h_val, sch)

    # CSV
    write_csv(wide, file.path(out_dir, paste0(stem, ".csv")))

    # LaTeX
    tex <- make_latex(wide, h_val, sch, n_targets)
    cat(as.character(tex), file = file.path(out_dir, paste0(stem, ".tex")))

    cat(sprintf("h=%2d, %-9s -> %s\n", h_val, sch, stem))
  }
}

cat("\nDone. Tables written to:", normalizePath(out_dir), "\n")
