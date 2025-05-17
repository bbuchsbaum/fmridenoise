context("ndx_rpca_temporal_components_multirun - Spike Mask")

.generate_spike_data <- function(T_run = 30, V = 20, N_runs = 2) {
  set.seed(123)
  base <- .generate_rpca_test_data(T_run = T_run, V = V, N_runs = N_runs)
  # Introduce large spikes in specific TRs
  spike_amplitude <- 50
  spike_tr_list <- list(c(10, 20), c(15))
  for (r in seq_len(N_runs)) {
    rows <- which(base$run_idx == r)
    spikes <- spike_tr_list[[r]]
    for (tr in spikes) {
      base$Y_residuals_cat[rows[tr], ] <- spike_amplitude
    }
  }
  base$expected_spikes_global <- c(10, 20, T_run + 15)
  base
}

test_that("spike_TR_mask flags injected spikes", {
  dat <- .generate_spike_data()
  user_opts <- list(
    rpca_merge_strategy = "concat_svd",
    rpca_max_iter = 5000,
    rpca_trace = FALSE,
    rpca_lambda_auto = FALSE,
    rpca_lambda_fixed = 0.1
  )

  res <- ndx_rpca_temporal_components_multirun(
    Y_residuals_cat = dat$Y_residuals_cat,
    run_idx = dat$run_idx,
    k_global_target = 2,
    user_options = user_opts
  )

  expect_true(is.list(res))
  expect_length(res$spike_TR_mask, dat$total_T)
  flagged <- which(res$spike_TR_mask)
  expect_true(all(dat$expected_spikes_global %in% flagged))
})
