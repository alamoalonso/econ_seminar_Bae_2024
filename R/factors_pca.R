#' PCA Factor Extraction
#'
#' Functions for extracting factors via Principal Component Analysis.
#'
#' @name factors_pca
NULL

#' Extract PCA factors
#'
#' @param X Matrix (T x N) of predictors, already standardized
#' @param k_max Maximum number of principal components to extract
#' @param config Configuration list
#' @return A list with components:
#'   - F: Matrix (T x k_max) of factor scores
#'   - loadings: Matrix (N x k_max) of factor loadings
#'   - pca_obj: The prcomp object
#' @export
#' @examples
#' X <- matrix(rnorm(100*50), 100, 50)
#' pca_result <- extract_pca(X, k_max = 5, config = config_us_default())
extract_pca <- function(X, k_max = 12, config = NULL) {
  # Input validation
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }

  if (nrow(X) < k_max) {
    stop(sprintf("Not enough observations (%d) to extract %d components", nrow(X), k_max))
  }

  if (ncol(X) < k_max) {
    log_warn(sprintf("Number of variables (%d) < k_max (%d). Using k_max = %d", ncol(X), k_max, ncol(X)), config)
    k_max <- ncol(X)
  }

  log_debug(sprintf("Extracting %d PCA factors from %d x %d matrix", k_max, nrow(X), ncol(X)), config)

  # Run PCA
  # X is already standardized, so center=TRUE, scale.=FALSE
  pca <- prcomp(X, center = isTRUE(config$pca_center), scale. = isTRUE(config$pca_scale))

  # Extract scores (T x k_max)
  F_hat <- pca$x[, 1:k_max, drop = FALSE]

  # Extract loadings (N x k_max)
  loadings <- pca$rotation[, 1:k_max, drop = FALSE]

  # Validate output dimensions
  if (nrow(F_hat) != nrow(X)) {
    stop(sprintf("PCA output dimension mismatch: expected %d rows, got %d", nrow(X), nrow(F_hat)))
  }
  if (ncol(F_hat) != k_max) {
    stop(sprintf("PCA output dimension mismatch: expected %d columns, got %d", k_max, ncol(F_hat)))
  }

  log_debug(sprintf("PCA complete: F is %d x %d", nrow(F_hat), ncol(F_hat)), config)

  list(
    F = F_hat,
    loadings = loadings,
    pca_obj = pca
  )
}
