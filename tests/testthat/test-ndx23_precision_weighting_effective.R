context("NDX-23 Precision Weighting Effectiveness")

# This test constructs a simple AR(2) signal with large spikes and checks
# that precision weights computed from the spike S matrix noticeably
# improve both AR coefficient estimation and regression betas.

test_that("precision weights greatly reduce spike influence", {
  set.seed(4242)
  n_tp <- 80
  phi_true <- c(0.6, -0.3)
  noise <- as.numeric(arima.sim(n = n_tp, model = list(ar = phi_true), sd = 0.1))

  X <- cbind(1, rnorm(n_tp))
  beta_true <- c(1, -1)

  Y <- as.vector(X %*% beta_true + noise)
  spike_idx <- c(25, 60)
  Y[spike_idx] <- Y[spike_idx] + c(20, -20)
  Y_mat <- matrix(Y, ncol = 1)

  S_mat <- matrix(0, nrow = n_tp, ncol = 1)
  S_mat[spike_idx, 1] <- c(20, -20)
  weights_mat <- ndx_precision_weights_from_S(S_mat)

  res_unw <- ndx_ar2_whitening(Y_mat, X, Y_mat, order = 2,
                               global_ar_on_design = FALSE,
                               verbose = FALSE)
  res_w <- ndx_ar2_whitening(Y_mat, X, Y_mat, order = 2,
                             global_ar_on_design = FALSE,
                             weights = weights_mat,
                             verbose = FALSE)

  err_unw <- sum(abs(res_unw$AR_coeffs_voxelwise - phi_true))
  err_w <- sum(abs(res_w$AR_coeffs_voxelwise - phi_true))
  expect_lt(err_w, err_unw)

  K_diag <- rep(0.05, ncol(X))
  betas_unw <- ndx_solve_anisotropic_ridge(res_unw$Y_whitened, X,
                                           K_penalty_diag = K_diag,
                                           na_mask = res_unw$na_mask)
  betas_w <- ndx_solve_anisotropic_ridge(res_w$Y_whitened, X,
                                         K_penalty_diag = K_diag,
                                         na_mask = res_w$na_mask,
                                         weights = weights_mat)
  diff_unw <- sum(abs(betas_unw[, 1] - beta_true))
  diff_w <- sum(abs(betas_w[, 1] - beta_true))
  expect_lt(diff_w, diff_unw)
})
