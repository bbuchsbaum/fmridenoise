context("ndx_estimate_initial_hrfs - Basic behavior")

set.seed(42)

# Simulate small fMRI dataset
n_time <- 40
n_voxels <- 4
TR_val <- 1
run_idx <- rep(1L, n_time)

y_signal <- sin(seq(0, 2*pi, length.out = n_time))
Y_fmri <- sapply(1:n_voxels, function(i) y_signal + rnorm(n_time, sd = 0.2))
Y_fmri <- as.matrix(Y_fmri)

# pass0 residuals with lower variance so some voxels have high R2
pass0_res <- matrix(rnorm(n_time * n_voxels, sd = 0.05), n_time, n_voxels)

# Two conditions with a few events each
events_df <- data.frame(
  onsets = c(5, 15, 25, 10, 20, 30),
  durations = rep(1, 6),
  condition = factor(rep(c("A", "B"), each = 3)),
  blockids = rep(1L, 6)
)

user_opts <- list(
  hrf_fir_taps = 6L,
  good_voxel_R2_threshold = 0.2,
  hrf_min_good_voxels = 1L,
  hrf_min_events_for_fir = 1L,
  hrf_low_event_threshold = 1L,
  lambda1_grid = 0.01,
  lambda2_grid = 0.01,
  cv_folds = 2L,
  hrf_cluster_method = "none",
  num_hrf_clusters = 1L
)

est_tbl <- ndx_estimate_initial_hrfs(Y_fmri, pass0_res, events_df,
                                     run_idx, TR_val,
                                     spike_TR_mask = NULL,
                                     user_options = user_opts)

test_that("returned tibble has one row per condition", {
  expect_s3_class(est_tbl, "tbl")
  expect_equal(nrow(est_tbl), length(unique(events_df$condition)))
})

test_that("hrf_estimate vectors have expected length", {
  expect_true(all(vapply(est_tbl$hrf_estimate, length, integer(1)) ==
                    user_opts$hrf_fir_taps))
})

test_that("num_effective_clusters attribute is 1 when clustering disabled", {
  expect_equal(attr(est_tbl, "num_effective_clusters"), 1L)
})


# High threshold: no good voxels remain -> NULL
user_opts_high_thr <- user_opts
user_opts_high_thr$good_voxel_R2_threshold <- 0.999

est_null <- ndx_estimate_initial_hrfs(Y_fmri, pass0_res, events_df,
                                      run_idx, TR_val,
                                      spike_TR_mask = NULL,
                                      user_options = user_opts_high_thr)

test_that("function returns NULL when no good voxels remain", {
  expect_null(est_null)
})
