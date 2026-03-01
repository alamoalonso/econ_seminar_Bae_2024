#!/usr/bin/env Rscript

# ============================================================================
# Table: MCS Results Summary  (Section 4.2.2)
# ============================================================================
# Produces the main MCS results tables from the paper outline.
#
# Input:
#   MCS_CSV  -- path to mcs_full_results.csv  (edit below)
#
# Output (outputs/tables/):
#   mcs_summary_<M0>_alpha<a>.{csv,tex}
#   One pair per M0_set_id x alpha combination.
#
# Table layout (12 rows x 6 columns):
#   Rows : all combinations of
#            scheme      in {recursive, rolling}
#            horizon     in {1, 12}
#            equation    in {DI, DIAR, DIAR-LAG}
#   Columns:
#     IR_{1-PLS}               inclusion rate of 1-PLS in M*
#     UR_{1-PLS}               uniqueness rate: share with M* = {1-PLS}
#     IR_{AR}                  inclusion rate of AR(BIC) in M*
#     E[|M*|]                  mean size of M*
#     E[|M*| | 1-PLS in M*]   conditional mean size given 1-PLS included
#     E[|M*| | AR in M*]       conditional mean size given AR included
#
# Main text table : alpha = 0.05
# Appendix tables : alpha = 0.10 and alpha = 0.01
# Separate tables : M0_set_id = "basic" (M_{0,1}) and "full" (M_{0,2})
# ============================================================================

if (basename(getwd()) == "scripts") setwd("..")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(knitr)
  library(kableExtra)
})

# ============================================================================
# Parameters  -- edit MCS_CSV if your run directory differs
# ============================================================================

MCS_CSV  <- "outputs/mcs_results_20260215_191157/mcs_full_results.csv"
HORIZONS <- c(1, 12)
SCHEMES  <- c("recursive", "rolling")
MODELS   <- c("DI", "DIAR", "DIAR-LAG")

# alpha for main text table; remaining alphas go to appendix tables
ALPHA_MAIN <- 0.05

out_dir <- file.path("outputs", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# Load and validate data
# ============================================================================

if (!file.exists(MCS_CSV)) {
  stop("MCS results file not found: ", MCS_CSV)
}

mcs_raw <- read_csv(MCS_CSV, show_col_types = FALSE)

required_cols <- c("series_id", "h", "scheme", "model_class",
                   "alpha", "M0_set_id", "n_methods_mcs",
                   "superior_set", "status")
missing_cols <- setdiff(required_cols, names(mcs_raw))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

# Keep only successfully completed MCS runs
mcs <- mcs_raw |>
  filter(status %in% c("ok", "fallback"),
         h %in% HORIZONS,
         scheme %in% SCHEMES,
         model_class %in% MODELS)

cat(sprintf("Loaded %d usable MCS rows  (%d total in file)\n",
            nrow(mcs), nrow(mcs_raw)))
cat(sprintf("  M0 sets : %s\n", paste(sort(unique(mcs$M0_set_id)), collapse = ", ")))
cat(sprintf("  alphas  : %s\n", paste(sort(unique(mcs$alpha)),     collapse = ", ")))
cat(sprintf("  horizons: %s\n", paste(sort(unique(mcs$h)),         collapse = ", ")))
cat(sprintf("  series  : %d\n\n", n_distinct(mcs$series_id)))

# ============================================================================
# Helper: check whether a method_id is in a semicolon-separated string
# ============================================================================

method_in_set <- function(method_id, set_str) {
  # set_str is a single character value like "recursive_PLS_grid_DI_k1;recursive_AR"
  # Returns TRUE if method_id appears as an exact token
  grepl(sprintf("(^|;)%s(;|$)", method_id), set_str)
}

# ============================================================================
# Core: compute summary statistics for one (M0_set_id, alpha) slice
# ============================================================================

compute_summary <- function(dat, m0, alp) {

  d <- dat |>
    filter(M0_set_id == m0, alpha == alp)

  if (nrow(d) == 0) return(NULL)

  # Build all 12 row conditions
  conditions <- expand_grid(
    scheme      = SCHEMES,
    h           = HORIZONS,
    model_class = MODELS
  )

  purrr::map_dfr(seq_len(nrow(conditions)), function(i) {

    sch <- conditions$scheme[i]
    h_  <- conditions$h[i]
    eq  <- conditions$model_class[i]

    rows <- d |>
      filter(scheme == sch, h == h_, model_class == eq)

    n_total <- nrow(rows)
    if (n_total == 0) {
      return(tibble(
        scheme = sch, h = h_, model_class = eq,
        IR_1PLS = NA_real_, UR_1PLS = NA_real_, IR_AR = NA_real_,
        mean_mcs_size = NA_real_,
        cond_size_1pls = NA_real_,
        cond_size_ar   = NA_real_,
        n_targets = 0L
      ))
    }

    # Method IDs for this slice
    id_1pls <- sprintf("%s_PLS_grid_%s_k1", sch, eq)
    id_ar   <- sprintf("%s_AR", sch)

    in_1pls <- method_in_set(id_1pls, rows$superior_set)
    in_ar   <- method_in_set(id_ar,   rows$superior_set)
    sz      <- rows$n_methods_mcs

    # Uniqueness: M* is exactly {1-PLS}  =>  in M* AND size == 1
    unique_1pls <- in_1pls & (sz == 1L)

    tibble(
      scheme      = sch,
      h           = h_,
      model_class = eq,
      # Inclusion rates
      IR_1PLS = mean(in_1pls),
      UR_1PLS = mean(unique_1pls),
      IR_AR   = mean(in_ar),
      # Mean set size (unconditional)
      mean_mcs_size  = mean(sz),
      # Conditional mean sizes
      cond_size_1pls = if (any(in_1pls))  mean(sz[in_1pls])  else NA_real_,
      cond_size_ar   = if (any(in_ar))    mean(sz[in_ar])    else NA_real_,
      n_targets = n_total
    )
  })
}

# ============================================================================
# Format summary tibble as a display-ready data frame
# ============================================================================

format_display <- function(summ) {
  summ |>
    mutate(
      scheme_label = if_else(scheme == "recursive", "Recursive", "Rolling"),
      h_label      = sprintf("$h = %d$", h),
      eq_label     = model_class,
      row_label    = sprintf("%s, %s", h_label, eq_label)
    ) |>
    arrange(
      factor(scheme, levels = SCHEMES),
      factor(h,           levels = HORIZONS),
      factor(model_class, levels = MODELS)
    ) |>
    transmute(
      scheme_label,
      `Horizon / Equation`             = row_label,
      `$IR_{1\\text{-PLS}}$`           = sprintf("%.2f", IR_1PLS),
      `$UR_{1\\text{-PLS}}$`           = sprintf("%.2f", UR_1PLS),
      `$IR_{AR}$`                      = sprintf("%.2f", IR_AR),
      `$\\mathbb{E}[|M^*|]$`           = sprintf("%.2f", mean_mcs_size),
      `$\\mathbb{E}[|M^*|\\mid 1\\text{-PLS}\\in M^*]$` =
        if_else(is.na(cond_size_1pls), "--", sprintf("%.2f", cond_size_1pls)),
      `$\\mathbb{E}[|M^*|\\mid AR\\in M^*]$` =
        if_else(is.na(cond_size_ar),   "--", sprintf("%.2f", cond_size_ar))
    )
}

# ============================================================================
# LaTeX table renderer
# ============================================================================

make_latex_mcs <- function(display_tbl, m0, alp, m0_size) {

  m0_label     <- if (m0 == "basic") "$M_{0,1}$" else "$M_{0,2}$"
  alpha_label  <- sprintf("$\\alpha = %.2f$", alp)
  scheme_sizes <- table(display_tbl$scheme_label)  # entries per scheme block

  # Drop the scheme_label grouping column before passing to kbl
  body <- display_tbl |> select(-scheme_label)

  col_names <- c("Horizon / Equation",
                 "$IR_{1\\text{-PLS}}$",
                 "$UR_{1\\text{-PLS}}$",
                 "$IR_{AR}$",
                 "$\\mathbb{E}[|M^*|]$",
                 "$\\mathbb{E}[|M^*|\\mid 1\\text{-PLS}\\in M^*]$",
                 "$\\mathbb{E}[|M^*|\\mid AR\\in M^*]$")

  caption <- sprintf(
    paste0(
      "MCS results summary. Candidate set: %s (%d methods). ",
      "Significance level: %s. ",
      "$IR$: inclusion rate in $M^*$ across targets. ",
      "$UR_{1\\text{-PLS}}$: share of targets with $M^* = \\{1\\text{-PLS}\\}$. ",
      "$\\mathbb{E}[|M^*|]$: mean size of the superior set. ",
      "Conditional mean sizes computed within the subset of targets ",
      "satisfying the stated condition."
    ),
    m0_label, m0_size, alpha_label
  )

  tbl <- kbl(
    body,
    format    = "latex",
    caption   = caption,
    col.names = col_names,
    booktabs  = TRUE,
    escape    = FALSE,
    align     = c("l", rep("r", 6))
  ) |>
    kable_styling(latex_options = c("hold_position", "scale_down")) |>
    add_header_above(
      c(" " = 1,
        "Inclusion / Uniqueness rates" = 3,
        "Mean size of $M^*$"           = 3),
      escape = FALSE
    )

  # Group rows by scheme
  schemes_in_order <- unique(display_tbl$scheme_label)
  row_counter <- 0L
  for (sch in schemes_in_order) {
    n <- sum(display_tbl$scheme_label == sch)
    tbl <- tbl |>
      group_rows(sch, row_counter + 1L, row_counter + n, bold = TRUE)
    row_counter <- row_counter + n
  }

  tbl
}

# ============================================================================
# Main loop: generate one table per (M0_set_id, alpha)
# ============================================================================

m0_sets   <- sort(unique(mcs$M0_set_id))
alphas    <- sort(unique(mcs$alpha))

# Candidate set sizes for caption
m0_sizes_tbl <- mcs |>
  group_by(M0_set_id) |>
  summarise(sz = max(n_methods_input, na.rm = TRUE), .groups = "drop")
m0_sizes <- setNames(m0_sizes_tbl$sz, m0_sizes_tbl$M0_set_id)

for (m0 in m0_sets) {
  cat(sprintf("=== M0 set: %s ===\n", m0))

  for (alp in alphas) {

    summ <- compute_summary(mcs, m0, alp)
    if (is.null(summ) || nrow(summ) == 0) {
      cat(sprintf("  alpha=%.2f : no data, skipped\n", alp))
      next
    }

    display <- format_display(summ)

    alpha_tag <- sub("\\.", "", sprintf("%.2f", alp))   # e.g. "005"
    stem      <- sprintf("mcs_summary_%s_alpha%s", m0, alpha_tag)

    # -- CSV --
    write_csv(summ, file.path(out_dir, paste0(stem, ".csv")))

    # -- LaTeX --
    m0_sz   <- m0_sizes[m0]
    tex     <- make_latex_mcs(display, m0, alp, m0_sz)
    tex_path <- file.path(out_dir, paste0(stem, ".tex"))
    cat(as.character(tex), file = tex_path)

    main_tag <- if (alp == ALPHA_MAIN) " [MAIN TEXT]" else " [APPENDIX]"
    cat(sprintf("  alpha=%.2f%s -> %s\n", alp, main_tag, stem))
  }
  cat("\n")
}

# ============================================================================
# Console preview: main text tables (alpha = 0.05)
# ============================================================================

cat(sprintf("=== Preview: alpha = %.2f ===\n\n", ALPHA_MAIN))

for (m0 in m0_sets) {
  summ <- compute_summary(mcs, m0, ALPHA_MAIN)
  if (is.null(summ)) next

  cat(sprintf("-- %s --\n", m0))
  print(
    summ |>
      select(scheme, h, model_class, IR_1PLS, UR_1PLS, IR_AR,
             mean_mcs_size, cond_size_1pls, cond_size_ar) |>
      mutate(across(where(is.numeric), ~ round(., 3))),
    n = Inf
  )
  cat("\n")
}

cat("Done. Tables written to:", normalizePath(out_dir), "\n")
