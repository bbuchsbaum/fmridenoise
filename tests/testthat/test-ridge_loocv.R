context("ndx_loocv_tune_lambda_ridge")

test_that("loocv lambda is selected from grid", {
  set.seed(42)
  n <- 30
  X <- cbind(1, rnorm(n))
  beta <- c(0.5, -1.2)
  Y <- X %*% beta + rnorm(n, sd = 0.05)
  lambda_grid <- c(0, 0.01, 0.1)
  lam <- ndx_loocv_tune_lambda_ridge(Y, X, lambda_grid)
  expect_true(lam %in% lambda_grid)
  expect_true(is.numeric(lam) && length(lam) == 1)
})
