#!/usr/bin/env Rscript
# Diagnostic Test: Check for Look-Ahead Bias in Factor Extraction
#
# This script traces through the forecasting pipeline to identify
# potential data leakage in factor estimation and regression fitting.

cat("=== LOOK-AHEAD BIAS DIAGNOSTIC ===\n\n")

# Source modules
source_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (f in source_files) source(f)

# Create minimal test case
cat("1. Creating minimal test dataset...\n")
set.seed(42)
T <- 20
N <- 10
dates <- seq.Date(as.Date("2010-01-01"), by = "month", length.out = T)

panel_final <- data.frame(
  date = dates,
  matrix(rnorm(T * N), T, N)
)
colnames(panel_final)[-1] <- paste0("X", 1:N)

# Create target
x_target <- rnorm(T)
panel_std1 <- data.frame(date = dates, TARGET = x_target)

cat(sprintf("   T = %d observations\n", T))
cat(sprintf("   Dates: %s to %s\n", dates[1], dates[T]))

# Construct h=1 target
cat("\n2. Constructing h=1 target...\n")
h <- 1
y_h <- construct_target_h(x_target, h)
cat("   y_h[t] = x[t+h]\n")
cat(sprintf("   y_h[1] = x[2] = %.3f\n", y_h[1]))
cat(sprintf("   y_h[2] = x[3] = %.3f\n", y_h[2]))
cat(sprintf("   y_h[T-1=%d] = x[T=%d] = %.3f\n", T-1, T, y_h[T-1]))
cat(sprintf("   y_h[T=%d] = NA (no x[T+1])\n", T))

targets_list <- list(
  TARGET = list(h1 = y_h)
)

# Simulate forecast at time T = 12
t_idx <- 12
t_origin <- dates[t_idx]

cat(sprintf("\n3. Simulating forecast at t_idx=%d (date=%s)...\n", t_idx, t_origin))
cat(sprintf("   Goal: forecast x[t+h] = x[%d]\n", t_idx + h))
cat(sprintf("   At time %d, we KNOW: x[1], ..., x[%d]\n", t_idx, t_idx))
cat(sprintf("   At time %d, we DO NOT KNOW: x[%d], ..., x[%d]\n", t_idx, t_idx+1, T))

# Current (buggy) implementation
cat("\n4. CURRENT implementation in extract_factors_at_origin:\n")
idx_win <- which(panel_final$date <= t_origin)
cat(sprintf("   idx_win <- which(date <= t_origin) = 1:%d\n", t_idx))

y_full <- targets_list[["TARGET"]][["h1"]]
y_win_CURRENT <- y_full[idx_win]

cat(sprintf("   y_win <- y_full[idx_win] = y_h[1:%d]\n", t_idx))
cat("   y_win contains:\n")
for (i in c(1, 2, t_idx-1, t_idx)) {
  if (i <= length(y_win_CURRENT)) {
    cat(sprintf("     y_h[%d] = x[%d] = %.3f", i, i+h, y_win_CURRENT[i]))
    if (i + h > t_idx) {
      cat(" <--- FUTURE VALUE! (x[%d] > x[%d])\n", i+h, t_idx)
    } else {
      cat(" (OK)\n")
    }
  }
}

cat(sprintf("\n   *** BUG: y_h[%d] = x[%d] is a FUTURE value! ***\n", t_idx, t_idx + h))
cat("   *** PLS uses this to estimate factor weights! ***\n")

# Correct implementation
cat("\n5. CORRECT implementation:\n")
valid_idx <- idx_win[idx_win + h <= t_idx]
cat(sprintf("   valid_idx <- idx_win[idx_win + h <= t_idx] = 1:%d\n", t_idx - h))

y_win_CORRECT <- y_full[valid_idx]
cat(sprintf("   y_win <- y_full[valid_idx] = y_h[1:%d]\n", length(valid_idx)))
cat("   y_win contains:\n")
for (i in seq_along(y_win_CORRECT)) {
  cat(sprintf("     y_h[%d] = x[%d] = %.3f (OK, x[%d] <= x[%d])\n",
              i, i+h, y_win_CORRECT[i], i+h, t_idx))
}

cat("\n6. Impact on PLS:\n")
cat("   CURRENT: PLS estimates weights using data up to y_h[t_idx] = x[t_idx+h]\n")
cat("            This includes the value we're trying to forecast!\n")
cat("   CORRECT: PLS should only use y_h[t] where t+h <= t_idx\n")
cat("            i.e., only use observed values at time t_idx\n")

# Check regression training indices
cat("\n7. Checking regression training indices:\n")
cat("   For DI model: y_h[idx] ~ F[idx]\n")

# Recursive window
idx_rec <- get_recursive_idx(t_idx)
cat(sprintf("   Recursive: idx_rec = 1:%d\n", t_idx))
cat(sprintf("   build_design_DI filters to: idx[!is.na(y_h[idx])]\n"))
cat(sprintf("   This includes indices up to t_idx=%d\n", t_idx))
cat(sprintf("   Training uses y_h[t_idx] = x[%d] (FUTURE!)\n", t_idx + h))

cat("\n8. Summary of Issues:\n")
cat("   Issue 1: PLS factor estimation uses future y values\n")
cat("            - At time T, y_win includes y_h[T] = x[T+h]\n")
cat("            - PLS weights are contaminated with look-ahead bias\n")
cat("\n")
cat("   Issue 2: Regression training may use future y values\n")
cat("            - DI/DIAR/DIAR-LAG fit using idx_rec or idx_roll\n")
cat("            - These indices include t_idx, giving y_h[t_idx] = x[t_idx+h]\n")
cat("            - Need to verify if filtering removes this\n")
cat("\n")

cat("\n=== DIAGNOSIS COMPLETE ===\n")
cat("\nCONCLUSION: LOOK-AHEAD BIAS DETECTED\n")
cat("  - PLS uses future y values when estimating factor weights\n")
cat("  - This explains monotonic improvement in k-PLS RMSE\n")
cat("  - Fix: Only use y_h[t] where t+h <= T when at forecast origin T\n")
