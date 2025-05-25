context("select_optimal_k_gdlite")

test_that("optimal K is selected on simple synthetic data", {
  set.seed(123)
  n_time <- 20
  run_idx <- rep(1:2, each = n_time/2)

  # base design: intercept and simple task regressor
  task <- rep(c(0, 1), length.out = n_time)
  X_base <- cbind(1, task)

  # dominant noise component (pc1)
  pc1 <- base::scale(seq_len(n_time), center = TRUE, scale = FALSE)
  pc1 <- pc1 / sqrt(sum(pc1^2))
  all_pcs <- cbind(pc1, rev(pc1))

  # generate data for 3 voxels
  betas <- matrix(c(1, 0.5, 0.2, 0.3, -0.1, 0.4), nrow = 2)
  Y_signal <- X_base %*% betas
  noise <- tcrossprod(pc1, c(2, 1, 1.5))
  Y_fmri <- Y_signal + noise

  good_mask <- rep(TRUE, 3)

  res <- select_optimal_k_gdlite(Y_fmri, X_base, all_pcs,
                                 run_idx, good_mask,
                                 K_max = 2, min_K = 0)

  expect_equal(res$K_star, 1)
  expect_true(is.matrix(res$selected_pcs))
  expect_equal(ncol(res$selected_pcs), 1)
  expect_equal(nrow(res$selected_pcs), n_time)
})
