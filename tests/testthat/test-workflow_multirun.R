context("NDX_Process_Subject - Multi-Run Workflow")

# This test ensures the workflow can handle multiple runs

test_that("Workflow runs with multiple runs and returns expected structure", {
  set.seed(99)
  TR_test <- 2.0
  n_time_per_run <- 20
  n_runs <- 2
  total_T <- n_time_per_run * n_runs
  n_vox <- 4

  Y_fmri <- matrix(rnorm(total_T * n_vox), nrow = total_T, ncol = n_vox)
  run_idx <- rep(as.integer(1:n_runs), each = n_time_per_run)
  motion_params <- matrix(rnorm(total_T * 3), nrow = total_T, ncol = 3)
  colnames(motion_params) <- paste0("mot", 1:3)

  events <- data.frame(
    onsets = c(5, 15, 5, 15) * TR_test,
    durations = rep(5, 4) * TR_test,
    condition = factor(rep(c("TaskA", "TaskB"), 2)),
    blockids = as.integer(rep(1:2, each = 2))
  )

  user_opts <- list(
    opts_pass0 = list(poly_degree = 1),
    opts_hrf = list(
      hrf_fir_taps = 6,
      hrf_fir_span_seconds = 12,
      good_voxel_R2_threshold = -Inf,
      lambda1_grid = c(0.1),
      lambda2_grid = c(0.1),
      cv_folds = 2,
      hrf_min_good_voxels = 1,
      hrf_cluster_method = "none",
      num_hrf_clusters = 1
    ),
    opts_rpca = list(
      k_global_target = 2,
      rpca_lambda_auto = FALSE,
      rpca_lambda_fixed = 0.1,
      rpca_merge_strategy = "concat_svd"
    ),
    opts_spectral = list(n_sine_candidates = 2, nyquist_guard_factor = 0.1),
    opts_whitening = list(global_ar_on_design = FALSE, max_ar_failures_prop = 0.5),
    opts_ridge = list(lambda_ridge = 0.5),
    max_passes = 2,
    min_des_gain_convergence = -Inf,
    min_rho_noise_projection_convergence = -Inf
  )

  res <- NDX_Process_Subject(
    Y_fmri = Y_fmri,
    events = events,
    motion_params = motion_params,
    run_idx = run_idx,
    TR = TR_test,
    user_options = user_opts,
    verbose = FALSE
  )

  expect_true(is.list(res))
  expect_equal(dim(res$Y_residuals_final_unwhitened), dim(Y_fmri))
  expect_equal(res$num_passes_completed, length(res$diagnostics_per_pass))
  expect_length(res$spike_TR_mask, total_T)
})
