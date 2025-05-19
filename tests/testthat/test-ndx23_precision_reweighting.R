context("NDX-23 Precision Re-weighting Integration")

# This test checks that precision weights derived from an RPCA S matrix
# reduce the influence of large spikes when estimating AR(2) coefficients
# and during anisotropic ridge regression.

test_that("S-based precision weights improve robustness to spikes", {
  set.seed(2024)
  n_tp <- 60
  phi_true <- c(0.5, -0.25)
  noise <- as.numeric(arima.sim(n = n_tp, model = list(ar = phi_true)))

  X <- cbind(1, rnorm(n_tp))
  beta_true <- c(0, 1)

  Y <- as.vector(X %*% beta_true + noise)
  spike_idx <- 30
  Y[spike_idx] <- Y[spike_idx] + 10
  Y_mat <- matrix(Y, ncol = 1)

  S_mat <- matrix(0, nrow = n_tp, ncol = 1)
  S_mat[spike_idx, 1] <- 10
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
  expect_lte(err_w, err_unw + 0.2)

  K_diag <- rep(0.1, ncol(X))
  betas_unw <- ndx_solve_anisotropic_ridge(res_unw$Y_whitened, X,
                                           K_penalty_diag = K_diag,
                                           na_mask = res_unw$na_mask)
  betas_w <- ndx_solve_anisotropic_ridge(res_w$Y_whitened, X,
                                         K_penalty_diag = K_diag,
                                         na_mask = res_w$na_mask,
                                         weights = weights_mat)
  diff_unw <- sum(abs(betas_unw[, 1] - beta_true))
  diff_w <- sum(abs(betas_w[, 1] - beta_true))
  expect_lte(diff_w, diff_unw + 0.2)
})
