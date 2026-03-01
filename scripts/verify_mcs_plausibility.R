# Verify MCS Results Plausibility
# This script checks if MCS results are consistent with RMSE rankings

library(dplyr)
library(tidyr)

# Load results
output_dir <- "outputs/20260215_164515"
rmse <- read.csv(file.path(output_dir, "rmse_results.csv"))
mcs <- read.csv(file.path(output_dir, "mcs_results.csv"))
forecasts <- read.csv(file.path(output_dir, "forecasts_long.csv"))

cat("=== MCS PLAUSIBILITY VERIFICATION ===\n\n")

# 1. Check that MCS status is "ok" for all tests
cat("1. MCS STATUS CHECK\n")
status_summary <- table(mcs$status)
print(status_summary)
cat("\n")

# 2. For each MCS test, check if the superior set contains the method with lowest MSE
cat("2. SUPERIOR SET vs RMSE RANKING\n\n")

# Create method_id in RMSE results to match MCS
rmse_with_method <- rmse %>%
  filter(!is.na(factor_method)) %>%
  mutate(
    method_id = case_when(
      k_mode == "grid" ~ paste(scheme, factor_spec_id, model, paste0("k", k), sep = "_"),
      k_mode == "dynamic" ~ paste(scheme, factor_spec_id, model, sep = "_"),
      TRUE ~ paste(scheme, "AR", sep = "_")
    )
  )

# Also add AR
ar_rmse <- rmse %>%
  filter(model == "AR") %>%
  mutate(method_id = paste(scheme, "AR", sep = "_"))

rmse_all <- bind_rows(rmse_with_method, ar_rmse)

# For each MCS row, check consistency
consistency_checks <- list()

for (i in 1:nrow(mcs)) {
  row <- mcs[i, ]

  # Parse included methods
  included <- strsplit(row$included_methods, ";")[[1]]
  superior <- strsplit(row$superior_set, ";")[[1]]

  # Get MSE for these methods from RMSE results
  mse_subset <- rmse_all %>%
    filter(
      series == row$series_id,
      h == row$h,
      scheme == row$scheme,
      model == row$model_class,
      method_id %in% included
    ) %>%
    select(method_id, mse, rmse_rel) %>%
    arrange(mse)

  if (nrow(mse_subset) == 0) next

  # Best method by MSE
  best_by_mse <- mse_subset$method_id[1]

  # Check if best method is in superior set
  is_consistent <- best_by_mse %in% superior

  consistency_checks[[length(consistency_checks) + 1]] <- data.frame(
    series_id = row$series_id,
    h = row$h,
    scheme = row$scheme,
    model_class = row$model_class,
    alpha = row$alpha,
    best_by_mse = best_by_mse,
    in_superior_set = is_consistent,
    superior_set = row$superior_set,
    n_mcs = row$n_methods_mcs,
    stringsAsFactors = FALSE
  )
}

consistency_df <- bind_rows(consistency_checks)

# Summary
cat("Consistency check: Best MSE method in MCS superior set?\n")
cat(sprintf("  Consistent: %d / %d (%.1f%%)\n",
    sum(consistency_df$in_superior_set),
    nrow(consistency_df),
    100 * mean(consistency_df$in_superior_set)))
cat("\n")

# Show inconsistent cases
inconsistent <- consistency_df %>% filter(!in_superior_set)
if (nrow(inconsistent) > 0) {
  cat("Inconsistent cases (best MSE not in MCS):\n")
  print(inconsistent)
} else {
  cat("All cases consistent - best MSE method always in MCS superior set.\n")
}
cat("\n")

# 3. Verify MCS p-values make sense
cat("3. P-VALUE ANALYSIS\n")
cat("Methods in superior set should have highest p-values.\n\n")

# Parse p-values for a sample case
sample_case <- mcs[1, ]
pvals_str <- sample_case$pvalues
pvals_pairs <- strsplit(pvals_str, ";")[[1]]
cat(sprintf("Sample case: %s, h=%d, %s, %s\n",
    sample_case$series_id, sample_case$h, sample_case$scheme, sample_case$model_class))
cat("P-values:\n")
for (p in pvals_pairs) {
  cat(sprintf("  %s\n", p))
}
cat("\n")

# 4. Check forecast loss matrix construction
cat("4. FORECAST LOSS MATRIX CHECK\n")
sample_forecasts <- forecasts %>%
  filter(series_id == "INDPRO", h == 1, scheme == "recursive", model_class == "DI") %>%
  group_by(method_id) %>%
  summarise(
    n_obs = n(),
    mean_se = mean((y_true - y_hat)^2, na.rm = TRUE),
    n_na = sum(is.na(y_hat)),
    .groups = "drop"
  )

cat("Forecasts for INDPRO, h=1, recursive, DI:\n")
print(sample_forecasts)
cat("\n")

# 5. Compare computed MSE with RMSE results
cat("5. MSE CONSISTENCY CHECK\n")
cat("Comparing MSE from forecasts vs rmse_results.csv:\n")

forecasts_mse <- forecasts %>%
  filter(series_id == "INDPRO", h == 1, scheme == "recursive") %>%
  group_by(method_id, model_class) %>%
  summarise(
    mse_from_forecasts = mean((y_true - y_hat)^2, na.rm = TRUE),
    .groups = "drop"
  )

rmse_mse <- rmse_all %>%
  filter(series == "INDPRO", h == 1, scheme == "recursive") %>%
  select(method_id, model, mse) %>%
  rename(model_class = model, mse_from_rmse = mse)

comparison <- left_join(forecasts_mse, rmse_mse, by = c("method_id", "model_class"))
comparison$diff <- abs(comparison$mse_from_forecasts - comparison$mse_from_rmse)
cat("MSE comparison (should be very small differences):\n")
print(comparison %>% select(method_id, mse_from_forecasts, mse_from_rmse, diff))
cat("\n")

# 6. Summary statistics
cat("=== SUMMARY ===\n")
cat(sprintf("Total MCS tests: %d\n", nrow(mcs)))
cat(sprintf("All status OK: %s\n", all(mcs$status == "ok")))
cat(sprintf("MCS size distribution:\n"))
print(table(mcs$n_methods_mcs))
cat("\n")

# Methods most frequently in MCS
superior_methods <- unlist(strsplit(mcs$superior_set, ";"))
method_freq <- table(superior_methods)
cat("Methods most frequently in MCS superior set:\n")
print(sort(method_freq, decreasing = TRUE)[1:min(10, length(method_freq))])

cat("\n=== VERIFICATION COMPLETE ===\n")
