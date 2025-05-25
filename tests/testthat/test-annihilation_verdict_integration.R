context("NDX_Process_Subject - Annihilation verdict integration")

test_that("Workflow outputs annihilation verdict when mode enabled", {
  set.seed(42)
  TR <- 1.5
  n_time <- 30
  n_vox <- 5
  Y <- matrix(rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)
  run_idx <- rep(1L, n_time)
  motion <- matrix(rnorm(n_time * 3), nrow = n_time, ncol = 3)
  colnames(motion) <- paste0("mot", 1:3)
  events <- data.frame(
    onsets = c(5, 15) * TR,
    durations = c(3, 3) * TR,
    condition = factor(c("A", "B")),
    blockids = 1L
  )

  opts <- list(
    opts_pass0 = list(poly_degree = 0),
    opts_hrf = list(
      hrf_fir_taps = 4,
      hrf_fir_span_seconds = 6,
      good_voxel_R2_threshold = -Inf,
      lambda1_grid = c(0.1),
      lambda2_grid = c(0.1),
      cv_folds = 2,
      hrf_min_good_voxels = 1,
      hrf_cluster_method = "none",
      num_hrf_clusters = 1
    ),
    opts_rpca = list(k_global_target = 2, rpca_lambda_auto = FALSE, rpca_lambda_fixed = 0.1, rpca_merge_strategy = "concat_svd"),
    opts_spectral = list(n_sine_candidates = 1, nyquist_guard_factor = 0.1),
    opts_whitening = list(global_ar_on_design = FALSE, max_ar_failures_prop = 0.5),
    opts_ridge = list(lambda_ridge = 0.5),
    opts_annihilation = list(
      annihilation_enable_mode = TRUE,
      annihilation_gdlite_k_max = 2,
      annihilation_gdlite_tsnr_thresh_noise_pool = 0.5,
      annihilation_gdlite_r2_thresh_noise_pool = 0.8,
      min_K_optimal_selection = 0
    ),
    max_passes = 1,
    min_des_gain_convergence = -Inf,
    min_rho_noise_projection_convergence = -Inf
  )

  res <- NDX_Process_Subject(
    Y_fmri = Y,
    events = events,
    motion_params = motion,
    run_idx = run_idx,
    TR = TR,
    user_options = opts,
    verbose = FALSE
  )

  expect_true(res$annihilation_mode_active)
  expect_true(!is.null(res$annihilation_verdict))
  expect_true(is.character(res$annihilation_verdict))
})
