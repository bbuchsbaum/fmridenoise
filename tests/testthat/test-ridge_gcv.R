context("anisotropic ridge GCV and variance")

test_that("ndx_estimate_res_var_whitened returns positive variance", {
  data <- .create_ridge_test_data(40, 2, 1, c(1, -1), noise_sd = 0.2)
  var_est <- ndx_estimate_res_var_whitened(data$Y, data$X)
  expect_true(is.numeric(var_est) && var_est > 0)
})

test_that("ndx_solve_anisotropic_ridge performs GCV tuning", {
  data <- .create_ridge_test_data(60, 2, 1, c(1, -1), noise_sd = 0.5)
  proj <- ndx_compute_projection_matrices(U_Noise = diag(2)[,1,drop=FALSE], n_regressors = 2)
  res_var <- ndx_estimate_res_var_whitened(data$Y, data$X)
  betas <- ndx_solve_anisotropic_ridge(data$Y, data$X,
    projection_mats = proj,
    lambda_values = list(lambda_parallel = 1, lambda_perp_signal = 0.1),
    gcv_lambda = TRUE,
    res_var_scale = res_var)$betas
  expect_equal(dim(betas), c(2,1))
})

test_that("warning when gcv_lambda TRUE but P_Noise missing", {
  data <- .create_ridge_test_data(40, 2, 1, c(1, -1), noise_sd = 0.5)
  proj <- ndx_compute_projection_matrices(U_Noise = NULL, n_regressors = 2)
  res_var <- ndx_estimate_res_var_whitened(data$Y, data$X)
  expect_warning(
    betas <- ndx_solve_anisotropic_ridge(data$Y, data$X,
      projection_mats = proj,
      lambda_values = list(lambda_parallel = 1, lambda_perp_signal = 0.1),
      gcv_lambda = TRUE,
      res_var_scale = res_var),
    "gcv_lambda = TRUE but projection_mats$P_Noise is missing", fixed = TRUE
  )
  expect_equal(dim(betas$betas), c(2,1))
})

test_that("ndx_gcv_tune_lambda_parallel handles NULL P_Noise", {
  data <- .create_ridge_test_data(30, 2, 1, c(1, -1), noise_sd = 0.2)
  lam <- ndx_gcv_tune_lambda_parallel(data$Y, data$X, NULL,
                                      lambda_grid = c(0.1, 1), lambda_ratio = 0.1)
  expect_true(is.numeric(lam) && length(lam) == 1)
  expect_true(lam %in% c(0.1, 1))
})
