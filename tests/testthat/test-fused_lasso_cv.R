context("cv_fusedlasso - cross validation works on synthetic data")

test_that("cv_fusedlasso returns best parameters from grid and finite errors", {
  set.seed(123)
  n <- 30
  p <- 6
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  beta_piecewise <- c(rep(1, 2), rep(0, 2), rep(-1, 2))
  y <- as.numeric(X %*% beta_piecewise + rnorm(n, sd = 0.1))

  block_ids <- rep(1:3, each = 10)
  lambda_grid <- 10^seq(-2, 0, length.out = 3)
  gamma_grid <- c(0.01, 0.1)
  k_folds <- 3

  res <- NULL
  expect_no_error({
    res <- cv_fusedlasso(y, X, lambda_grid, gamma_grid, k_folds, block_ids)
  })

  expect_true(res$best_lambda %in% lambda_grid,
              label = "best_lambda should come from lambda_grid")
  expect_true(res$best_gamma %in% gamma_grid,
              label = "best_gamma should come from gamma_grid")
  expect_equal(dim(res$cv_error_matrix), c(length(lambda_grid), length(gamma_grid)),
               label = "cv_error_matrix should have correct dimensions")
  expect_true(all(is.finite(res$cv_error_matrix)),
              label = "cv_error_matrix should have finite values")
})
