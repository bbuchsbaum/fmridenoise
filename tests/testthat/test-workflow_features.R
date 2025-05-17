context("NDX_Process_Subject - Additional Features")

# Test 1: spike_TR_mask propagation and RPCA spikes

test_that("NDX_Process_Subject propagates spike mask and detects spikes", {
  set.seed(123)
  dat <- .generate_spike_data(T_run = 30, V = 5, N_runs = 2)
  TR_test <- 2.0
  total_T <- dat$total_T

  Y_fmri_test <- dat$Y_residuals_cat
  run_idx_test <- dat$run_idx
  motion_params_test <- matrix(0, nrow = total_T, ncol = 3)
  colnames(motion_params_test) <- paste0("mot", 1:3)

  events_test <- data.frame(
    onsets = c(5, 15, 5, 15) * TR_test,
    durations = rep(2 * TR_test, 4),
    condition = factor(c("TaskA", "TaskB", "TaskA", "TaskB")),
    blockids = c(1, 1, 2, 2)
  )

  initial_mask <- rep(FALSE, total_T)
  initial_mask[1] <- TRUE

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
      rpca_lambda_fixed = 0.1
    ),
    opts_spectral = list(n_sine_candidates = 2, nyquist_guard_factor = 0.1),
    opts_whitening = list(global_ar_on_design = FALSE, max_ar_failures_prop = 0.5),
    opts_ridge = list(lambda_ridge = 0.5),
    max_passes = 1,
    min_des_gain_convergence = -Inf,
    min_rho_noise_projection_convergence = -Inf
  )

  res <- NDX_Process_Subject(
    Y_fmri = Y_fmri_test,
    events = events_test,
    motion_params = motion_params_test,
    run_idx = run_idx_test,
    TR = TR_test,
    spike_TR_mask = initial_mask,
    user_options = user_opts,
    verbose = FALSE
  )

  expect_length(res$spike_TR_mask, total_T)
  expect_type(res$spike_TR_mask, "logical")
  expect_true(all(which(initial_mask) %in% which(res$spike_TR_mask)))
  flagged <- which(res$spike_TR_mask)
  expect_true(all(dat$expected_spikes_global %in% flagged))
})

# Test 2: Annihilation mode activates GLMdenoise-Lite PCs

test_that("NDX_Process_Subject produces gdlite PCs in annihilation mode", {
  set.seed(42)
  TR_test <- 2.0
  n_time <- 40
  n_voxels <- 5
  Y_fmri <- matrix(rnorm(n_time * n_voxels), nrow = n_time, ncol = n_voxels)
  run_idx <- rep(1, n_time)
  motion_params <- matrix(0, nrow = n_time, ncol = 2)
  colnames(motion_params) <- paste0("mot", 1:2)

  events <- data.frame(
    onsets = c(10, 20),
    durations = c(5, 5),
    condition = factor(c("TaskA", "TaskB")),
    blockids = 1
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
    opts_rpca = list(k_global_target = 2, rpca_lambda_auto = FALSE, rpca_lambda_fixed = 0.1),
    opts_spectral = list(n_sine_candidates = 2, nyquist_guard_factor = 0.1),
    opts_whitening = list(global_ar_on_design = FALSE, max_ar_failures_prop = 0.5),
    opts_ridge = list(lambda_ridge = 0.5),
    opts_annihilation = list(
      annihilation_enable_mode = TRUE,
      annihilation_gdlite_k_max = 2,
      annihilation_gdlite_tsnr_thresh_noise_pool = 0,
      min_K_optimal_selection = 1
    ),
    max_passes = 1,
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

  expect_true(res$annihilation_mode_active)
  expect_true(is.matrix(res$gdlite_pcs))
  expect_true(ncol(res$gdlite_pcs) > 0)
  expect_true(any(grepl("^gdlite_pc_", colnames(res$X_full_design_final))))
})

