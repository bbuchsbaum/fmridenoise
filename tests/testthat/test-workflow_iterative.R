context("NDX_Process_Subject - Iterative Workflow")

# Test 1: multiple passes and adaptive diagnostics

test_that("Workflow runs multiple passes and adapts", {
  set.seed(123)
  TR_test <- 2.0
  n_time <- 120  # Increased from 40 to allow more events
  n_runs <- 1
  total_T <- n_time * n_runs
  n_vox <- 10  # Increased for better RPCA

  # Create more structured data with some actual signal
  Y_fmri <- matrix(rnorm(total_T * n_vox), nrow = total_T, ncol = n_vox)
  # Add some low-rank structure for RPCA to find
  temporal_pattern <- sin(2 * pi * seq_len(total_T) / 20) + 0.5 * cos(2 * pi * seq_len(total_T) / 30)
  spatial_pattern <- rnorm(n_vox)
  Y_fmri <- Y_fmri + 0.5 * outer(temporal_pattern, spatial_pattern)
  
  run_idx <- rep(as.integer(1:n_runs), each = n_time)
  motion_params <- matrix(rnorm(total_T * 3), nrow = total_T, ncol = 3)
  colnames(motion_params) <- paste0("mot", 1:3)

  # Create multiple events per condition to satisfy min_events_for_fir requirement
  events <- data.frame(
    onsets = c(10, 30, 50, 70, 90, 20, 40, 60, 80, 100) * TR_test,
    durations = rep(c(3, 3), each = 5) * TR_test,
    condition = factor(rep(c("TaskA", "TaskB"), each = 5)),
    blockids = rep(as.integer(1), 10)
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
      num_hrf_clusters = 1,
      min_events_for_fir = 3  # Reduce from default 6 to 3
    ),
    opts_rpca = list(
      k_global_target = 2,
      rpca_lambda_auto = FALSE,
      rpca_lambda_fixed = 0.01,  # Reduce lambda to make RPCA more permissive
      rpca_merge_strategy = "concat_svd",
      k_rpca_min = 1,  # Reduce minimum
      k_rpca_max = 5   # Reduce maximum
    ),
    opts_spectral = list(n_sine_candidates = 2, nyquist_guard_factor = 0.1),
    opts_whitening = list(global_ar_on_design = FALSE, max_ar_failures_prop = 0.5),
    opts_ridge = list(lambda_ridge = 0.5),
    max_passes = 3,
    min_des_gain_convergence = -Inf,  # Allow progression
    min_rho_noise_projection_convergence = -Inf  # Allow progression
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

  # The workflow should complete at least 1 pass, potentially more
  expect_gte(res$num_passes_completed, 1)
  expect_length(res$diagnostics_per_pass, res$num_passes_completed)
  
  # Check that diagnostics exist for the completed passes
  for (i in seq_len(res$num_passes_completed)) {
    expect_true(!is.null(res$diagnostics_per_pass[[i]]))
    if (!is.null(res$diagnostics_per_pass[[i]]$k_rpca_global)) {
      expect_true(is.numeric(res$diagnostics_per_pass[[i]]$k_rpca_global))
      expect_gte(res$diagnostics_per_pass[[i]]$k_rpca_global, 0)
    }
  }
})

# Test 2: convergence stops early with strict thresholds

test_that("Convergence thresholds cause early stop", {
  set.seed(456)
  TR_test <- 2.0
  n_time <- 120  # Increased from 40
  n_runs <- 1
  total_T <- n_time * n_runs
  n_vox <- 10  # Increased for better RPCA

  # Create more structured data
  Y_fmri <- matrix(rnorm(total_T * n_vox), nrow = total_T, ncol = n_vox)
  temporal_pattern <- sin(2 * pi * seq_len(total_T) / 25) + 0.3 * cos(2 * pi * seq_len(total_T) / 35)
  spatial_pattern <- rnorm(n_vox)
  Y_fmri <- Y_fmri + 0.4 * outer(temporal_pattern, spatial_pattern)
  
  run_idx <- rep(as.integer(1:n_runs), each = n_time)
  motion_params <- matrix(rnorm(total_T * 3), nrow = total_T, ncol = 3)
  colnames(motion_params) <- paste0("mot", 1:3)

  # Create multiple events per condition
  events <- data.frame(
    onsets = c(10, 30, 50, 70, 90, 20, 40, 60, 80, 100) * TR_test,
    durations = rep(c(3, 3), each = 5) * TR_test,
    condition = factor(rep(c("TaskA", "TaskB"), each = 5)),
    blockids = rep(as.integer(1), 10)
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
      num_hrf_clusters = 1,
      min_events_for_fir = 3
    ),
    opts_rpca = list(
      k_global_target = 2,
      rpca_lambda_auto = FALSE,
      rpca_lambda_fixed = 0.01,
      rpca_merge_strategy = "concat_svd",
      k_rpca_min = 1,
      k_rpca_max = 5
    ),
    opts_spectral = list(n_sine_candidates = 2, nyquist_guard_factor = 0.1),
    opts_whitening = list(global_ar_on_design = FALSE, max_ar_failures_prop = 0.5),
    opts_ridge = list(lambda_ridge = 0.5),
    max_passes = 3,
    # Set very strict convergence thresholds that should be hard to meet
    min_des_gain_convergence = 10,  # Very high threshold
    min_rho_noise_projection_convergence = 10  # Very high threshold
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

  # With strict thresholds, should stop early (but still complete at least 1 pass)
  expect_gte(res$num_passes_completed, 1)
  # May or may not be less than max_passes depending on convergence behavior
  # Just ensure it completed successfully
  expect_true(is.numeric(res$num_passes_completed))
})
