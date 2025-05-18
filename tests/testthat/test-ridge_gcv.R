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
    res_var_scale = res_var)
  expect_equal(dim(betas), c(2,1))
})
