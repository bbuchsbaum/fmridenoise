context("HRF estimation helper functions")

library(matrixStats)

# Test validate_hrf_inputs


test_that("validate_hrf_inputs merges defaults and validates data", {
  set.seed(1)
  Y <- matrix(rnorm(40), nrow = 20, ncol = 2)
  pass0 <- matrix(rnorm(40), nrow = 20, ncol = 2)
  events_df <- data.frame(
    onsets = c(0, 10),
    durations = c(2, 2),
    condition = c("A", "B"),
    blockids = c(1, 1)
  )
  run_idx <- rep(1, 20)
  user_opts <- list(hrf_fir_taps = 8, good_voxel_R2_threshold = 0.2)

  val <- validate_hrf_inputs(Y, pass0, events_df, run_idx, TR = 2,
                             spike_TR_mask = NULL, user_opts)

  expect_true(is.list(val))
  expect_equal(val$n_timepoints, 20)
  expect_equal(length(val$validated_spike_TR_mask), 20)
  expect_true(is.integer(val$user_options$hrf_fir_taps))
  expect_equal(val$user_options$hrf_fir_taps, 8L)
  expect_equal(val$user_options$good_voxel_R2_threshold, 0.2)

  expect_error(validate_hrf_inputs(Y, pass0[1:10,], events_df, run_idx, 2, NULL, user_opts),
               "pass0_residuals must be a numeric matrix")
  expect_error(validate_hrf_inputs(Y, pass0, events_df, run_idx[1:10], 2, NULL, user_opts),
               "run_idx must be a numeric vector")
})


test_that("prepare_hrf_response_data selects voxels and respects spike mask", {
  set.seed(42)
  Y <- cbind(seq_len(10), seq_len(10) + 1)
  pass0 <- matrix(rnorm(20, sd = 0.01), nrow = 10, ncol = 2)
  run_idx <- rep(1, 10)
  spike_mask <- c(TRUE, TRUE, rep(FALSE, 8))
  user_opts <- list(good_voxel_R2_threshold = 0.1, hrf_min_good_voxels = 1)

  res <- prepare_hrf_response_data(Y, pass0, run_idx, spike_mask, user_opts)

  expect_true(is.list(res))
  expect_equal(res$n_good_voxels, 2)
  expect_equal(length(res$ybar_clean), 8)
  expect_equal(res$block_ids_for_cv, run_idx[!spike_mask])
  expect_equal(res$good_voxels_indices, c(1, 2))
})


test_that("prepare_hrf_response_data returns NULL with few good voxels", {
  set.seed(99)
  Y <- matrix(rnorm(20), nrow = 10, ncol = 2)
  pass0 <- matrix(rnorm(20, sd = 5), nrow = 10, ncol = 2)
  run_idx <- rep(1, 10)
  spike_mask <- rep(FALSE, 10)
  user_opts <- list(good_voxel_R2_threshold = 0.9, hrf_min_good_voxels = 3)

  expect_warning(res_null <- prepare_hrf_response_data(Y, pass0, run_idx, spike_mask, user_opts),
                 "Insufficient good voxels")
  expect_null(res_null)
})


