#' @keywords internal
# NO @export for cv_fusedlasso
cv_fusedlasso <- function(y, X, lambda_grid, gamma_grid, k_folds, block_ids) {
  # Minimal block-wise cross-validation for fusedlasso
  if (length(lambda_grid) == 0) lambda_grid <- 0
  if (length(gamma_grid) == 0) gamma_grid <- 0

  unique_blocks <- unique(block_ids)
  k_eff <- min(length(unique_blocks), max(1L, k_folds))
  fold_ids <- rep(seq_len(k_eff), length.out = length(unique_blocks))
  names(fold_ids) <- unique_blocks

  build_D <- function(p) {
    if (p > 1) {
      diff(diag(p), differences = 1)
    } else {
      matrix(0, nrow = 0, ncol = 1)
    }
  }

  D <- build_D(ncol(X))

  best_lambda <- lambda_grid[1]
  best_gamma <- gamma_grid[1]
  best_score <- Inf

  for (lam in lambda_grid) {
    for (gam in gamma_grid) {
      mse_folds <- c()
      for (fold in seq_len(k_eff)) {
        test_blocks <- unique_blocks[fold_ids == fold]
        train_idx <- !(block_ids %in% test_blocks)
        test_idx <- (block_ids %in% test_blocks)
        if (sum(train_idx) == 0 || sum(test_idx) == 0) next

        fit <- tryCatch(
          genlasso::fusedlasso(y[train_idx], X[train_idx, , drop = FALSE], D = D, gamma = gam),
          error = function(e) NULL
        )
        if (is.null(fit)) next

        beta <- tryCatch(stats::coef(fit, lambda = lam)$beta, error = function(e) NULL)
        if (is.null(beta)) next

        pred <- as.numeric(X[test_idx, , drop = FALSE] %*% beta)
        mse_folds <- c(mse_folds, mean((y[test_idx] - pred)^2))
      }

      if (length(mse_folds) == 0) {
        # Fallback: use training error on full data
        fit <- tryCatch(
          genlasso::fusedlasso(y, X, D = D, gamma = gam),
          error = function(e) NULL
        )
        if (!is.null(fit)) {
          beta <- tryCatch(stats::coef(fit, lambda = lam)$beta, error = function(e) NULL)
          if (!is.null(beta)) {
            pred <- as.numeric(X %*% beta)
            mse_folds <- mean((y - pred)^2)
          }
        }
      } else {
        mse_folds <- mean(mse_folds)
      }

      if (length(mse_folds) > 0 && is.finite(mse_folds) && mse_folds < best_score) {
        best_score <- mse_folds
        best_lambda <- lam
        best_gamma <- gam
      }
    }
  }

  list(best_lambda = best_lambda, best_gamma = best_gamma)
}
