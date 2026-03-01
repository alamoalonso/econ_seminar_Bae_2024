# Econometrics Seminar Paper — Working Outline & Reporting Plan (Persistent Notes)

This file captures the agreed paper structure and the current plan for reporting and interpreting results, especially for the Model Confidence Set (MCS) extension. It is intended as a stable reference while you write the paper.

Last updated: 2026-02-23


## 1) Paper structure (LaTeX section plan)

You will implement the following structure in LaTeX (no abstract).

### \section{Introduction}

Single section without subheadings. Paragraph flow:

1) Motivation for factor models: **N large vs T** forecasting setting; why dimension reduction is useful.  
2) Selective literature positioning: most relevant theoretical + empirical + practical application strands.  
3) Smooth transition to **Bae (2024)**: her research question, brief methodology, and headline finding (neutral phrasing).  
4) Transition to **this paper’s** research question: whether there is **statistical evidence** that **1-PLS** is dominant among alternatives; explain what you add (replication/interpretation + MCS extension).  
5) Roadmap of remaining sections.


### \section{Theory}

#### \subsection{Approximate dynamic factor model}
- State the general approximate (dynamic) factor model.
- State the general forecasting equation (as in the presentation), without introducing DI/DIAR/DIAR-LAG yet.

#### \subsection{Factor estimation methods}

##### \subsubsection{Principal Component Analysis (PCA)}
- Precise definition of PCA factor estimation used in the paper.

##### \subsubsection{Partial Least Squares (PLS)}
- Precise definition of the PLS factor estimation procedure consistent with the presentation notation.
- Emphasize: PLS yields **target-specific** supervised factors; PCA yields **unsupervised** factors.

#### \subsection{Parameter choice for $k$}

Planned structure:

- **2.3:** k-choice problem (why $k$ matters; bias–variance / signal–noise motivation; why decision rules exist). Discuss theoretical properties of decision rules (like consistency with regards to the true number of latent factors) with reference to papers investigating these properties.
- **2.3.1:** Fixed-k methods
- **2.3.2:** Bai & Ng for PCA (BN-BIC)
- **2.3.3:** Onatski for PLS (ON)

Key commitment: decision rules will be defined **mathematically precisely**.

#### \subsection{Model Confidence Set}
- Formal and mathematically precise introduction of the MCS procedure:
  - Candidate set $M_0$, loss and loss differentials, EPA null, test statistic used, bootstrap details, elimination rule, output $M^*$ at confidence level $1-\alpha$.


### \section{Data and Empirical Design}

#### \subsection{Dataset}
- Document dataset used, all relevant transformations, sample design, horizons, and recursive vs rolling scheme definition. State explicitly that, from the 148 forecasting target used in Bae's empirical study, we only use X (to be specified later) and refer to an overview of targets used in the Appendix.

#### \subsection{Forecasting equations}
- Formally introduce **DI**, **DIAR**, **DIAR-LAG** as special cases of the general forecasting equation by specifying admissible lag structures.
- Define **BIC** mathematically (general definition) and document how it is used to select lag orders (e.g., $p^*$ and $m^*$).
- Briefly motivate how/why the forecasting equations differ and why DI may be less competitive.

#### \subsection{Evaluation}
- Define the forecast loss (MSE / RMSE) precisely.
- Define the **AR(BIC)** benchmark model precisely (including BIC-based order selection).
- Define **Relative MSE/RMSE** against AR(BIC) (as used in the replication).


### \section{Empirical Results}

#### \subsection{Replication}
- State explicitly you focus on a subset of Bae’s empirical evidence.

##### \subsubsection{Mean RMSE across $k$}
- Replicate the “mean RMSE across $k$” evidence (Figures 1–2 style argument).
- Extend with heterogeneity evidence:
  - Individual-target plots (already prepared).
  - A table counting how often each $k\in\{1,\dots,12\}$ attains minimum OOS RMSE, with:
    - Rows: equation × horizon, structured as:
      - DI: h=1, h=6, h=12
      - DIAR: …
      - DIAR-LAG: …
    - show one scheme in main text; other scheme in appendix.

##### \subsubsection{Diebold–Mariano tests against Benchmark}
- Replicate DM tests against AR(BIC).
- Explain the key limitation: DM here is benchmark-relative and does not establish “best among alternatives,” motivating the MCS extension.

#### \subsection{Extension: MCS}

##### \subsubsection{Experimental setup}
Planned design:
- Run MCS **separately** for:
  - target variables (all),
  - horizons **h ∈ {1, 12}**,
  - schemes **{recursive, rolling}**,
  - forecasting equations **{DI, DIAR, DIAR-LAG}**.
- Candidate sets:
  - $M_{0,1} = \{\text{AR(BIC)},\ 1\text{-PLS},\ 1\text{-PCA},\ \text{PLS(ON)},\ \text{PCA(BN-BIC)}\}$
  - $M_{0,2} = \{\text{AR(BIC)},\ 1\text{-PLS},\ 2\text{-PLS},\ 3\text{-PLS},\ 1\text{-PCA},\ 2\text{-PCA},\ 3\text{-PCA},\ \text{PLS(ON)},\ \text{PCA(BN-BIC)}\}$
- Rationale:
  - Include ON (for PLS) and BN-BIC (for PCA) as the best-performing non-fixed decision rules in Bae (2024).
  - Add fixed-k alternatives (k=2,3) in $M_{0,2}$ since fixed-k methods frequently attain minimum OOS RMSE in the replication.
  - Include AR(BIC) (benchmark) in both $M_0$ sets to track how often factor methods collectively outperform the benchmark.

##### \subsubsection{MCS Results}
Main reporting strategy (to avoid overloading the reader):
- Produce one **main results table** for each $M_0$ (i.e., one for $M_{0,1}$ and one for $M_{0,2}$).
- Fix the reported significance level in the main text to **α = 0.05**.
- Put tables for **α ∈ {0.10, 0.01}** into the appendix.

**Rows (12 conditions):**
- All combinations of:
  - scheme ∈ {recursive, rolling}
  - horizon ∈ {1, 12}
  - forecasting equation ∈ {DI, DIAR, DIAR-LAG}

**Columns (fixed measures):**
1) $IR_{1\text{-PLS}}$  
2) $UR_{1\text{-PLS}}$  
3) $IR_{AR}$  
4) $\mathbb{E}[|M^*|]$  
5) $\mathbb{E}[|M^*|\mid 1\text{-PLS} \in M^*]$  
6) $\mathbb{E}[|M^*|\mid AR \in M^*]$

Definitions:
- $IR_{1\text{-PLS}}$: Inclusion rate of 1-PLS in $M^*$ across targets within a given condition.  
- $UR_{1\text{-PLS}}$: Uniqueness rate of 1-PLS, defined as share of targets with $M^* = \{1\text{-PLS}\}$.  
- $IR_{AR}$: Inclusion rate of AR(BIC) in $M^*$.  
- $\mathbb{E}[|M^*|]$: Mean size of the superior set.  
- Conditional mean set sizes computed within the subset of targets where the condition holds.

Appendix policy:
- Any additional interesting summaries (plots, pairwise survival frequencies, target-category breakdowns, etc.) will go to the appendix only.

##### \subsubsection{MCS Interpretation}
To be written once the computed tables are analyzed. Planned interpretation logic (neutral, observational):

- Interpret $IR_{1\text{-PLS}}$ as “frequency with which 1-PLS remains statistically indistinguishable from the best within the candidate set.”  
- Interpret $UR_{1\text{-PLS}}$ as a stronger indicator of “uniquely best” performance (singleton superior set).  
- Interpret $IR_{AR}$ as diagnostic of whether the benchmark remains competitive and whether factor methods provide statistically detectable improvements under the MCS procedure.  
- Interpret $\mathbb{E}[|M^*|]$ (and conditional mean sizes) as indicators of how often the data support a small set of clearly superior models versus a large set of statistically indistinguishable models.

Concrete wording and claims will be finalized after inspecting the computed tables (α=0.05 main text; α=0.10/0.01 appendix).


### \section{Conclusion}
- Summarize replication findings and what they do (and do not) support regarding “1-PLS dominance.”
- Summarize what the MCS extension adds relative to benchmark-only testing.
- Provide a brief critical assessment and limitations (e.g., dependence/bootstrapping choices, candidate-set sensitivity via $M_{0,1}$ vs $M_{0,2}$), keeping tone neutral and research-like.


## 2) Implementation notes / guardrails

- Tone: observational, research-like; avoid judgmental phrasing.
- Main text prioritizes: formal definitions + key empirical evidence + compact MCS summaries.
- Appendices: alternative α tables, secondary plots, extra granular breakdowns.
