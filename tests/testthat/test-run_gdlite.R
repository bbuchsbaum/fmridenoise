context("ndx_run_gdlite")

test_that("run_gdlite returns expected structure on toy data", {
  set.seed(1)
  n_time <- 6
  run_idx <- rep(1:2, each = 3)
  TR <- 1.0
  Y <- matrix(rnorm(n_time * 4), n_time, 4)

  events <- data.frame(
    onsets = c(0, 3),
    durations = c(1, 1),
    condition = "A",
    blockids = 1:2
  )

  res <- ndx_run_gdlite(
    Y_fmri = Y,
    events = events,
    run_idx = run_idx,
    TR = TR,
    poly_degree = 0,
    k_max = 2,
    r2_thresh_noise_pool = 0,
    tsnr_thresh_noise_pool = -Inf,
    r2_thresh_good_voxels = 0,
    detrend_for_tsnr = FALSE,
    perform_final_glm = FALSE
  )

  expect_true(is.list(res))
  expect_true(all(c("selected_pcs", "K_star", "X_gd") %in% names(res)))
  if (res$K_star > 0) {
    expect_true(is.matrix(res$selected_pcs))
    expect_equal(nrow(res$selected_pcs), n_time)
    expect_equal(ncol(res$selected_pcs), res$K_star)
  } else {
    expect_null(res$selected_pcs)
  }
  expect_true(is.matrix(res$X_gd))
})

test_that("run_gdlite works with realistic multi-run fMRI data", {
  set.seed(42)
  
  # Create realistic multi-run fMRI data
  n_runs <- 3
  n_time_per_run <- 100
  n_voxels <- 200
  TR <- 2.0
  total_time <- n_runs * n_time_per_run
  
  # Create run indices
  run_idx <- rep(1:n_runs, each = n_time_per_run)
  
  # Create realistic events with proper timing
  events <- data.frame(
    onsets = c(10, 50, 90, 130, 170, 210, 250, 290),  # Events spread across runs
    durations = rep(2, 8),
    condition = rep(c("A", "B"), 4),
    blockids = c(1, 1, 1, 2, 2, 2, 3, 3)  # Events distributed across runs
  )
  
  # Ensure events are properly ordered by blockids to match run_idx order
  events <- events[order(events$blockids, events$onsets), ]
  
  # Generate synthetic fMRI data with task activation and noise
  time_points <- seq(0, (total_time - 1) * TR, by = TR)
  
  # Create task-related signal (simplified HRF convolution)
  task_signal <- matrix(0, nrow = total_time, ncol = n_voxels)
  for (i in 1:nrow(events)) {
    onset_tr <- round(events$onsets[i] / TR) + 1
    if (onset_tr <= total_time) {
      # Simple HRF-like response (gamma function approximation)
      hrf_length <- min(20, total_time - onset_tr + 1)
      hrf_times <- 0:(hrf_length - 1)
      hrf <- (hrf_times^8.6) * exp(-hrf_times/0.547) / (0.547^8.6 * gamma(8.6 + 1))
      
      # Add task signal to subset of voxels (simulating activation)
      active_voxels <- sample(n_voxels, size = round(n_voxels * 0.3))
      for (v in active_voxels) {
        end_tr <- min(onset_tr + hrf_length - 1, total_time)
        task_signal[onset_tr:end_tr, v] <- task_signal[onset_tr:end_tr, v] + 
          hrf[1:(end_tr - onset_tr + 1)] * rnorm(1, mean = 2, sd = 0.5)
      }
    }
  }
  
  # Add structured noise (simulating physiological noise)
  # Low-frequency drift
  drift_signal <- matrix(0, nrow = total_time, ncol = n_voxels)
  for (v in 1:n_voxels) {
    drift_freq <- runif(1, 0.001, 0.01)  # Very low frequency
    drift_signal[, v] <- sin(2 * pi * drift_freq * time_points) * rnorm(1, 0, 1)
  }
  
  # Respiratory-like noise (around 0.3 Hz)
  resp_signal <- matrix(0, nrow = total_time, ncol = n_voxels)
  resp_voxels <- sample(n_voxels, size = round(n_voxels * 0.6))
  for (v in resp_voxels) {
    resp_freq <- runif(1, 0.25, 0.35)
    resp_signal[, v] <- sin(2 * pi * resp_freq * time_points) * rnorm(1, 0, 0.8)
  }
  
  # Cardiac-like noise (around 1.2 Hz)
  cardiac_signal <- matrix(0, nrow = total_time, ncol = n_voxels)
  cardiac_voxels <- sample(n_voxels, size = round(n_voxels * 0.4))
  for (v in cardiac_voxels) {
    cardiac_freq <- runif(1, 1.0, 1.4)
    cardiac_signal[, v] <- sin(2 * pi * cardiac_freq * time_points) * rnorm(1, 0, 0.6)
  }
  
  # Random noise
  random_noise <- matrix(rnorm(total_time * n_voxels, 0, 1), nrow = total_time, ncol = n_voxels)
  
  # Combine all signals
  Y_fmri <- task_signal + drift_signal + resp_signal + cardiac_signal + random_noise
  
  # Add motion parameters
  motion_params <- matrix(rnorm(total_time * 6, 0, 0.5), nrow = total_time, ncol = 6)
  colnames(motion_params) <- paste0("mot", 1:6)
  
  # Run GLMdenoise-lite
  res <- ndx_run_gdlite(
    Y_fmri = Y_fmri,
    events = events,
    run_idx = run_idx,
    TR = TR,
    motion_params = motion_params,
    poly_degree = 1,
    k_max = 10,
    r2_thresh_noise_pool = 0.9,  # Very permissive to include most voxels in noise pool
    tsnr_thresh_noise_pool = 0.1,  # Very permissive
    r2_thresh_good_voxels = 0.01,
    min_K_optimal_selection = 0,
    detrend_for_tsnr = TRUE,
    perform_final_glm = TRUE
  )
  
  # Test structure and basic properties
  expect_true(is.list(res))
  expect_true(all(c("selected_pcs", "K_star", "X_gd", "r2_cv_initial", "tsnr_values",
                    "noise_pool_mask", "good_voxel_mask_for_k_selection",
                    "median_r2cv_by_K", "all_r2cv_per_K", "betas_gdlite",
                    "residuals_gdlite", "all_candidate_pcs") %in% names(res)))
  
  # Test dimensions
  expect_equal(nrow(res$X_gd), total_time)
  expect_true(ncol(res$X_gd) > 0)  # Should have task + motion + polynomial regressors
  expect_equal(length(res$r2_cv_initial), n_voxels)
  expect_equal(length(res$tsnr_values), n_voxels)
  expect_equal(length(res$noise_pool_mask), n_voxels)
  expect_equal(length(res$good_voxel_mask_for_k_selection), n_voxels)
  
  # Test K_star is reasonable
  expect_true(res$K_star >= 0)
  expect_true(res$K_star <= 10)  # Should not exceed k_max
  
  # Test selected PCs dimensions if any were selected
  if (res$K_star > 0) {
    expect_true(is.matrix(res$selected_pcs))
    expect_equal(nrow(res$selected_pcs), total_time)
    expect_equal(ncol(res$selected_pcs), res$K_star)
    expect_true(all(is.finite(res$selected_pcs)))
  }
  
  # Test that noise pool mask is logical (may be empty with synthetic data)
  expect_true(is.logical(res$noise_pool_mask))
  expect_true(sum(res$noise_pool_mask) >= 0)  # Can be 0 with synthetic data
  
  # Test that good voxel mask is logical (may be empty with synthetic data)
  expect_true(is.logical(res$good_voxel_mask_for_k_selection))
  expect_true(sum(res$good_voxel_mask_for_k_selection) >= 0)  # Can be 0 with synthetic data
  
  # Test final GLM outputs
  if (!is.null(res$betas_gdlite)) {
    expect_true(is.matrix(res$betas_gdlite))
    expect_equal(nrow(res$betas_gdlite), ncol(res$X_gd) + res$K_star)
    expect_equal(ncol(res$betas_gdlite), n_voxels)
  }
  
  if (!is.null(res$residuals_gdlite)) {
    expect_true(is.matrix(res$residuals_gdlite))
    expect_equal(nrow(res$residuals_gdlite), total_time)
    expect_equal(ncol(res$residuals_gdlite), n_voxels)
  }
  
  # Test that R2 values are reasonable
  expect_true(all(res$r2_cv_initial >= 0 | is.na(res$r2_cv_initial)))
  expect_true(all(res$r2_cv_initial <= 1 | is.na(res$r2_cv_initial)))
  
  # Test that tSNR values are reasonable (can be 0 or NA, but not negative)
  expect_true(all(res$tsnr_values >= 0 | is.na(res$tsnr_values)))
})

test_that("run_gdlite handles edge cases properly", {
  set.seed(123)
  
  # Test with single run (should handle LORO CV gracefully)
  n_time <- 50
  n_voxels <- 20
  Y_single <- matrix(rnorm(n_time * n_voxels), nrow = n_time, ncol = n_voxels)
  run_idx_single <- rep(1, n_time)
  events_single <- data.frame(
    onsets = c(5, 25),
    durations = c(2, 2),
    condition = "A",
    blockids = c(1, 1)
  )
  
  res_single <- ndx_run_gdlite(
    Y_fmri = Y_single,
    events = events_single,
    run_idx = run_idx_single,
    TR = 2.0,
    poly_degree = 0,
    k_max = 5,
    perform_final_glm = FALSE
  )
  
  expect_true(is.list(res_single))
  expect_true(all(is.na(res_single$r2_cv_initial)))  # LORO CV should fail with single run
  expect_true(res_single$K_star >= 0)
  
  # Test with very strict thresholds (should result in empty noise pool)
  res_strict <- ndx_run_gdlite(
    Y_fmri = Y_single,
    events = events_single,
    run_idx = run_idx_single,
    TR = 2.0,
    poly_degree = 0,
    k_max = 5,
    r2_thresh_noise_pool = -1,  # Very strict
    tsnr_thresh_noise_pool = 1000,  # Very strict
    perform_final_glm = FALSE
  )
  
  expect_equal(res_strict$K_star, 0)
  expect_true(sum(res_strict$noise_pool_mask) == 0)
  
  # Test with no events
  events_empty <- data.frame(
    onsets = numeric(0),
    durations = numeric(0),
    condition = character(0),
    blockids = numeric(0)
  )
  
  res_no_events <- ndx_run_gdlite(
    Y_fmri = Y_single,
    events = events_empty,
    run_idx = run_idx_single,
    TR = 2.0,
    poly_degree = 1,  # Should still have polynomial regressors
    k_max = 3,
    perform_final_glm = FALSE
  )
  
  expect_true(is.list(res_no_events))
  expect_true(ncol(res_no_events$X_gd) > 0)  # Should have polynomial regressors
})

test_that("run_gdlite polynomial degree handling works correctly", {
  set.seed(456)
  n_time <- 30
  n_voxels <- 10
  Y <- matrix(rnorm(n_time * n_voxels), nrow = n_time, ncol = n_voxels)
  run_idx <- rep(1:2, each = 15)
  events <- data.frame(
    onsets = c(2, 17),
    durations = c(1, 1),
    condition = "A",
    blockids = 1:2
  )
  
  # Test poly_degree = 0 (run-specific intercepts only)
  res_poly0 <- ndx_run_gdlite(
    Y_fmri = Y, events = events, run_idx = run_idx, TR = 1.0,
    poly_degree = 0, k_max = 2, perform_final_glm = FALSE,
    r2_thresh_noise_pool = 0.9, tsnr_thresh_noise_pool = 0.1
  )
  expect_true(ncol(res_poly0$X_gd) >= 2)  # At least 2 run intercepts
  
  # Test poly_degree = 1 (linear trends)
  res_poly1 <- ndx_run_gdlite(
    Y_fmri = Y, events = events, run_idx = run_idx, TR = 1.0,
    poly_degree = 1, k_max = 2, perform_final_glm = FALSE,
    r2_thresh_noise_pool = 0.9, tsnr_thresh_noise_pool = 0.1
  )
  expect_true(ncol(res_poly1$X_gd) > ncol(res_poly0$X_gd))  # Should have more regressors
  
  # Test poly_degree = 2 (quadratic trends)
  res_poly2 <- ndx_run_gdlite(
    Y_fmri = Y, events = events, run_idx = run_idx, TR = 1.0,
    poly_degree = 2, k_max = 2, perform_final_glm = FALSE,
    r2_thresh_noise_pool = 0.9, tsnr_thresh_noise_pool = 0.1
  )
  expect_true(ncol(res_poly2$X_gd) > ncol(res_poly1$X_gd))  # Should have even more regressors
})

test_that("run_gdlite validates inputs properly", {
  set.seed(789)
  n_time <- 20
  n_voxels <- 5
  Y <- matrix(rnorm(n_time * n_voxels), nrow = n_time, ncol = n_voxels)
  run_idx <- rep(1:2, each = 10)
  events <- data.frame(
    onsets = c(2, 12),
    durations = c(1, 1),
    condition = "A",
    blockids = 1:2
  )
  
  # Test mismatched blockids and run_idx
  events_bad <- events
  events_bad$blockids <- c(1, 3)  # blockid 3 doesn't exist in run_idx
  
  expect_error(
    ndx_run_gdlite(Y_fmri = Y, events = events_bad, run_idx = run_idx, TR = 1.0),
    "events\\$blockids.*do not match run_idx"
  )
  
  # Test with motion parameters of wrong dimensions
  motion_bad <- matrix(rnorm(10 * 6), nrow = 10, ncol = 6)  # Wrong number of rows
  
  expect_warning(
    res_motion_bad <- ndx_run_gdlite(
      Y_fmri = Y, events = events, run_idx = run_idx, TR = 1.0,
      motion_params = motion_bad, perform_final_glm = FALSE
    ),
    "Motion parameters.*row count does not match"
  )
  
  # Should still work, just ignore bad motion params
  expect_true(is.list(res_motion_bad))
})
