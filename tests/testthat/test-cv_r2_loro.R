context("cv_r2_loro")

test_that("perfectly predicted data yields R2 of one", {
  n_runs <- 2
  n_per_run <- 4
  run_idx <- rep(1:n_runs, each = n_per_run)
  time <- 1:(n_runs * n_per_run)

  X <- cbind(1, sin(time))
  betas <- matrix(c(2, 3), nrow = 2)
  Y <- X %*% betas

  r2 <- cv_r2_loro(Y, X, run_idx)
  expect_equal(r2, rep(1, ncol(Y)))
})

test_that("design matrix with zero columns returns zeros", {
  Y <- matrix(rnorm(10 * 3), nrow = 10)
  X <- matrix(nrow = 10, ncol = 0)
  run_idx <- rep(1:2, each = 5)
  expect_warning(r2 <- cv_r2_loro(Y, X, run_idx))
  expect_equal(r2, rep(0, ncol(Y)))
})

test_that("voxels with zero variance give R2 of zero", {
  run_idx <- rep(1:2, each = 3)
  X <- matrix(1, nrow = 6, ncol = 1)
  Y <- cbind(rep(5, 6), rnorm(6))
  r2 <- cv_r2_loro(Y, X, run_idx)
  expect_equal(r2[1], 0)
  expect_true(r2[2] >= 0 && r2[2] <= 1)
})
