You are working inside my econometrics forecasting replication repository (R project). Do NOT modify any files yet. First, assess the status quo and produce a detailed implementation plan that I can verify before you start coding.

Goal: extend the forecasting pipeline so it supports (1) PCA factors with Bai–Ng BN–BIC (BIC3) decision rule for selecting the number of factors k, and (2) PLS factors with Onatski (2010) “ON” decision rule for selecting k. We also must reliably save the produced forecasts for each method/horizon/target/time to outputs so they can be used as inputs to a Model Confidence Set (MCS) procedure later.

Constraints / expectations:
- Step 1: Inspect the repository structure, identify where factor estimation, k-selection, forecast regression, evaluation, and output-writing currently live.
- Step 2: Identify what “methods” exist now (e.g., PCA, PLS, 1-PLS, existing decision rules like BN-p1/p2/p3, BIC, GP, etc.), how these are configured (config file? function arguments? method registry?), and how forecasts are currently generated/stored.
- Step 3: Propose a clean architecture to add decision rules as pluggable components, and to make method definitions like PCA(BN-BIC) and PLS(ON) first-class, reproducible entries.
- Step 4: Provide a plan to implement (later) with minimal refactoring but high reliability: where to insert k-selection logic, how to compute k recursively at each forecast origin, and how to store outputs in a stable format for MCS.
- Step 5: Identify potential ambiguities (data standardization, rolling/recursive windows, indexing with horizon h, maximum k, numerical stability) and propose default choices aligned with Bae (2024) style.
- Deliverable NOW: a written plan + file/function touch-points, including proposed function signatures, required new modules/files (if any), output schema, and a checklist for verification.

Statistical details you MUST use (do not guess from memory; implement exactly as specified):

A) PCA with Bai–Ng BN–BIC decision rule (BIC3 from Bai & Ng 2002)
- At each forecast origin t (recursive/expanding window), determine k by minimizing BN’s BIC3 over k in {0,1,…,k_max}.
- Use the PCA reconstruction residual variance:
  Let X be the N x T_t training predictor matrix at origin t (after standardization used elsewhere in repo).
  For each k:
    - Estimate k factors via PCA on X.
    - Compute residuals E_k = X - \hatΛ_k \hatF_k' (the rank-k approximation implied by PCA).
    - Define V(k) = (1/(N*T_t)) * sum_{i=1..N} sum_{τ=1..T_t} E_{k,iτ}^2.
  Also define \hatσ^2 = V(k_max) (common practice in Bai–Ng IC definitions; confirm if repo already does something else; note it in plan).
  Then compute:
    BIC3(k) = ln(V(k)) + k * \hatσ^2 * ((N + T_t - k)/(N*T_t)) * ln(N*T_t)
  (If the repo’s conventions use V(k) without log, document it, but for BN-BIC we want the BIC3 form above. If you find the repo already has BN ICs implemented, map exactly to their internal form and document equivalence.)
  Choose k_hat = argmin_k BIC3(k).

B) PLS with Onatski (2010) “ON” decision rule
- ON is based purely on eigenvalues of the predictor covariance, not y, so it avoids the “V(k) for PLS” ambiguity.
- At each forecast origin t, compute k_hat_ON(t) from the eigen-spectrum of the sample covariance of X in the training window:
  Let X be N x T_t training predictor matrix (same X used to extract factors).
  Compute sample covariance in the “variables dimension”:
    Σ_hat = (1/T_t) * X X' (N x N).
  Let λ_1 >= λ_2 >= ... >= λ_N be eigenvalues of Σ_hat.
  Choose r_max such that 2*r_max + 1 <= min(N, T_t). Use a repo-level default like r_max = 12 unless the repo already defines r_max; in the plan, propose a robust rule and make it configurable.
  Define:
    w = 2^(2/3) / (2^(2/3) - 1)
    u_hat = w * λ_{r_max+1} + (1 - w) * λ_{2*r_max+1}
  Define δ(t) = max(N^(-2/5), T_t^(-2/5)) (Onatski’s suggested shrinking threshold for finite samples; make it configurable).
  Then:
    k_hat_ON = #{ i : λ_i > (1 + δ(t)) * u_hat }
- Use k_hat_ON as the number of PLS components/factors to extract (run the repo’s iterative PLS algorithm up to k_hat_ON and then use these factors in the forecasting regression).

C) Recursive updating / where to compute k
- Both BN-BIC and ON must be applied within the recursive forecasting loop: at each origin t, compute k_hat(t) on the training sample available at t, then estimate factors using that k_hat(t), then estimate the forecast regression and generate the forecast for horizon h.

Outputs / saving forecasts for MCS
- We need a deterministic forecast output table per run. Propose an output schema (CSV/parquet/RDS) with at least:
  method_id (e.g., “PCA_BN-BIC”, “PLS_ON”, “1-PLS”), target_id, horizon_h, forecast_origin_t (date/index), y_true (realized), y_hat, training_window_info (start/end index), k_hat (selected factors), and any other metadata required to reproduce results (seed, standardization regime, lag order p if applicable).
- Ensure that forecasts from all methods are stored in the SAME standardized schema so MCS can read them without special casing.
- Plan should include where in the pipeline these outputs are written, how filenames are constructed, and how to avoid accidental overwrites (use a run_id / timestamp / git hash if available).

Repository assessment instructions
- Start by listing key folders and scripts; identify entrypoints (e.g., run_forecasts.R, main.R, pipeline.R, etc.).
- Find the code that currently:
  1) constructs X and y for each target/horizon
  2) loops over forecast origins
  3) estimates factors (PCA/PLS)
  4) selects k (if any rules exist)
  5) estimates the forecasting regression
  6) writes forecast outputs
- Summarize current behavior, gaps relative to required behavior, and propose minimal changes.

Deliverable: write a structured plan with:
1) Status quo summary (what exists, where)
2) Proposed design (modules/functions)
3) Step-by-step implementation plan (ordered tasks)
4) Proposed output schema + example row
5) Verification checklist (unit tests / sanity checks)
6) Any risks/edge cases + mitigations

Do NOT implement code yet.
