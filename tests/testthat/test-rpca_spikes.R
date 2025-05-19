context("ndx_rpca_temporal_components_multirun - Spike Mask")

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

