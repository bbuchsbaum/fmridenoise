context("anisotropic ridge - additional tests")

# helper to build simple dataset
.make_simple_data <- function(n_tp, X, betas_true, noise_sd=0.1) {
  Y <- X %*% betas_true + matrix(rnorm(n_tp * ncol(betas_true), sd = noise_sd), n_tp)
  return(list(Y=Y, X=X))
}


test_that("projection-based penalties equal explicit penalties", {
  set.seed(42)
  n_tp <- 30
  X <- matrix(rnorm(n_tp*4), n_tp, 4)
  betas_true <- matrix(c(1, -1, 2, -2), nrow=4, ncol=1)
  data <- .make_simple_data(n_tp, X, betas_true, noise_sd=0.1)

  U_noise <- diag(4)[,1:2]
  proj <- ndx_compute_projection_matrices(U_Noise = U_noise, n_regressors = 4)
  lambda_vals <- list(lambda_parallel = 3, lambda_perp_signal = 0.5)

  betas_proj <- ndx_solve_anisotropic_ridge(data$Y, data$X,
                                           projection_mats = proj,
                                           lambda_values = lambda_vals)

  K_diag <- c(rep(lambda_vals$lambda_parallel,2),
              rep(lambda_vals$lambda_perp_signal,2))
  betas_diag <- ndx_solve_anisotropic_ridge(data$Y, data$X,
                                            K_penalty_diag = K_diag)

  dimnames(betas_proj) <- NULL
  dimnames(betas_diag) <- NULL
  expect_equal(betas_proj, betas_diag, tolerance = 1e-6)
})


test_that("anisotropic penalties preserve signal while shrinking noise", {
  set.seed(99)
  n_tp <- 40
  X <- matrix(rnorm(n_tp*2), n_tp, 2)
  betas_true <- matrix(c(0, 3), nrow=2, ncol=1)
  data <- .make_simple_data(n_tp, X, betas_true, noise_sd=0.2)

  U_noise <- diag(2)[,1, drop=FALSE]
  proj <- ndx_compute_projection_matrices(U_Noise = U_noise, n_regressors = 2)
  lambda_vals <- list(lambda_parallel = 100, lambda_perp_signal = 0.1)

  betas_aniso <- ndx_solve_anisotropic_ridge(data$Y, data$X,
                                            projection_mats = proj,
                                            lambda_values = lambda_vals)
  K_diag_iso <- rep(lambda_vals$lambda_parallel, 2)
  betas_iso <- ndx_solve_anisotropic_ridge(data$Y, data$X,
                                           K_penalty_diag = K_diag_iso)

  error_aniso <- abs(betas_aniso[2,1] - betas_true[2,1])
  error_iso <- abs(betas_iso[2,1] - betas_true[2,1])

  expect_lt(abs(betas_aniso[1,1]), 0.2)
  expect_lt(error_aniso, error_iso)
})

