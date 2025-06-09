context("ndx_run_sprint1 - Workflow Integration Tests")

# Load helper functions
source("helper-rpca_data.R")

# Generate realistic test data using proven helper function
TR_test <- 2.0
n_time_per_run_test <- 60  # Longer runs for more realistic data
n_runs_test <- 2           # Multiple runs
n_voxels_test <- 30        # More voxels for better component detection

# Use the helper function that generates data with known low-rank structure
realistic_data <- .generate_rpca_test_data(
  T_run = n_time_per_run_test,
  V = n_voxels_test, 
  N_runs = n_runs_test
)

Y_fmri_test <- realistic_data$Y_residuals_cat
run_idx_test <- realistic_data$run_idx
total_timepoints_test <- realistic_data$total_T

# Create realistic motion parameters
motion_params_test <- matrix(rnorm(total_timepoints_test * 6, sd = 0.1), 
                            nrow = total_timepoints_test, ncol = 6)
colnames(motion_params_test) <- paste0("mot", 1:6)

# Create events with sufficient repetitions for HRF estimation
events_per_run <- 6  # More events for better HRF estimation
onsets_per_run <- seq(10, n_time_per_run_test - 10, length.out = events_per_run) * TR_test

events_list <- list()
for (run in 1:n_runs_test) {
  run_events <- data.frame(
    onsets = onsets_per_run,
    durations = rep(4 * TR_test, events_per_run),  # 4 TR duration
    condition = factor(rep(c("TaskA", "TaskB"), length.out = events_per_run)),
    blockids = as.integer(rep(run, events_per_run))
  )
  events_list[[run]] <- run_events
}
events_test <- do.call(rbind, events_list)

# User options optimized for realistic data
user_options_test <- list(
  opts_pass0 = list(
    poly_degree = 2  # Higher degree for longer runs
  ),
  opts_hrf = list(
    hrf_fir_taps = 8,
    hrf_fir_span_seconds = 16,
    good_voxel_R2_threshold = 0.01,  # Reasonable threshold
    lambda1_grid = c(0.01, 0.1, 1.0),  # Broader grid
    lambda2_grid = c(0.01, 0.1, 1.0),
    cv_folds = 3,
    hrf_min_good_voxels = 5,
    hrf_cluster_method = "none",
    num_hrf_clusters = 1,
    hrf_min_events_for_fir = 3  # Lower threshold since we have 6 events per condition
  ),
  opts_rpca = list(
    k_global_target = 4,  # Higher k for more complex data
    rpca_lambda_auto = FALSE,
    rpca_lambda_fixed = 0.001  # Very low lambda for better detection
  ),
  opts_spectral = list(
    n_sine_candidates = 5,  # More candidates
    nyquist_guard_factor = 0.8,
    k_tapers = 3,
    nw = 2
  ),
  opts_whitening = list(
    global_ar_on_design = TRUE,  # More realistic setting
    max_ar_failures_prop = 0.3
  ),
  opts_ridge = list(
    lambda_ridge = 0.1,  # Lower regularization
    anisotropic_ridge_enable = TRUE  # Test anisotropic ridge
  ),
  task_regressor_names_for_extraction = c("task_TaskA", "task_TaskB"),
  max_passes = 3,  # More passes to test iterative behavior
  min_des_gain_convergence = 0.001,  # Realistic convergence criteria
  min_rho_noise_projection_convergence = 0.01
)

test_that("NDX_Process_Subject completes workflow successfully with realistic data", {
  
  workflow_output <- NULL
  expect_no_error({
    workflow_output <- NDX_Process_Subject(
      Y_fmri = Y_fmri_test,
      events = events_test,
      motion_params = motion_params_test,
      run_idx = run_idx_test,
      TR = TR_test,
      user_options = user_options_test,
      verbose = FALSE
    )
  })
  
  expect_true(is.list(workflow_output), "Workflow output should be a list")
  
  expected_names <- c(
    "annihilation_mode_active", "ar_coeffs_voxelwise", "beta_history_per_pass", 
    "diagnostics_per_pass", "estimated_hrfs", "final_task_betas", "ljung_box_p", 
    "na_mask_whitening", "num_hrf_clusters", "num_passes_completed", "pass0_residuals", 
    "pass0_vars", "rpca_components", "S_matrix_rpca_final", "spectral_sines", 
    "spike_TR_mask", "X_full_design_final", "Y_residuals_final_unwhitened"
  )
  expect_named(workflow_output, expected = expected_names, ignore.order = TRUE)
  
  # Basic dimension checks for key matrix outputs (these should always be present)
  expect_equal(dim(workflow_output$Y_residuals_final_unwhitened), dim(Y_fmri_test))
  expect_true(is.numeric(workflow_output$pass0_vars) && length(workflow_output$pass0_vars) == 1)
  expect_equal(dim(workflow_output$pass0_residuals), dim(Y_fmri_test))
  
  # Workflow should complete at least one pass
  expect_true(workflow_output$num_passes_completed >= 1)
  
  # Spike mask should be properly initialized
  expect_length(workflow_output$spike_TR_mask, total_timepoints_test)
  expect_type(workflow_output$spike_TR_mask, "logical")
  
  # Check diagnostics_per_pass structure
  expect_true(is.list(workflow_output$diagnostics_per_pass))
  expect_length(workflow_output$diagnostics_per_pass, workflow_output$num_passes_completed)
  
  if (workflow_output$num_passes_completed > 0) {
    expect_true(all(sapply(workflow_output$diagnostics_per_pass, function(p) "DES" %in% names(p))))
    
    # DES should be calculated (may be NA if workflow fails early, but should be present)
    first_pass_diag <- workflow_output$diagnostics_per_pass[[1]]
    expect_true("DES" %in% names(first_pass_diag), "DES should be present in diagnostics")
  }
  
  expect_true(is.list(workflow_output$beta_history_per_pass))
  expect_length(workflow_output$beta_history_per_pass, workflow_output$num_passes_completed)
  
  # Check AR coefficients if present
  if (!is.null(workflow_output$ar_coeffs_voxelwise)) {
    expect_true(is.matrix(workflow_output$ar_coeffs_voxelwise))
    expect_equal(ncol(workflow_output$ar_coeffs_voxelwise), n_voxels_test)
    expect_equal(nrow(workflow_output$ar_coeffs_voxelwise), 2)  # AR(2) coefficients
  }

  # Check RPCA spike detection if present
  if (!is.null(workflow_output$S_matrix_rpca_final)) {
    expect_equal(dim(workflow_output$S_matrix_rpca_final), dim(Y_fmri_test))
  }
  
  # Annihilation mode should be inactive by default
  expect_false(workflow_output$annihilation_mode_active)
})

test_that("NDX_Process_Subject produces meaningful outputs when components succeed", {
  # This test focuses on checking that when components do work, they produce sensible outputs
  # We'll use a simpler configuration that's more likely to succeed
  
  simple_options <- list(
    opts_pass0 = list(poly_degree = 1),
    opts_hrf = list(
      hrf_fir_taps = 6, hrf_fir_span_seconds = 12, good_voxel_R2_threshold = -Inf,
      lambda1_grid = c(0.1), lambda2_grid = c(0.1), cv_folds = 2,
      hrf_min_good_voxels = 1, hrf_cluster_method = "none", num_hrf_clusters = 1,
      hrf_min_events_for_fir = 2
    ),
    opts_rpca = list(k_global_target = 2, rpca_lambda_auto = FALSE, rpca_lambda_fixed = 0.001),
    opts_spectral = list(n_sine_candidates = 3, nyquist_guard_factor = 0.8),
    opts_whitening = list(global_ar_on_design = FALSE, max_ar_failures_prop = 0.5),
    opts_ridge = list(lambda_ridge = 0.1, anisotropic_ridge_enable = FALSE),  # Use simpler isotropic ridge
    task_regressor_names_for_extraction = c("task_TaskA", "task_TaskB"),
    max_passes = 2, min_des_gain_convergence = -Inf, min_rho_noise_projection_convergence = -Inf
  )
  
  result <- NDX_Process_Subject(
    Y_fmri = Y_fmri_test,
    events = events_test,
    motion_params = motion_params_test,
    run_idx = run_idx_test,
    TR = TR_test,
    user_options = simple_options,
    verbose = FALSE
  )
  
  # Check that the workflow produces valid outputs
  expect_true(result$num_passes_completed >= 1)
  
  # If RPCA succeeds, check its outputs
  if (!is.null(result$rpca_components) && ncol(result$rpca_components) > 0) {
    expect_equal(nrow(result$rpca_components), total_timepoints_test)
    expect_true(ncol(result$rpca_components) <= user_options_test$opts_rpca$k_global_target)
  }
  
  # If spectral analysis succeeds, check its outputs  
  if (!is.null(result$spectral_sines) && ncol(result$spectral_sines) > 0) {
    expect_equal(nrow(result$spectral_sines), total_timepoints_test)
    expect_true(ncol(result$spectral_sines) %% 2 == 0)  # Should be pairs of sin/cos
  }
  
  # If design matrix is created, check its structure
  if (!is.null(result$X_full_design_final)) {
    expect_equal(nrow(result$X_full_design_final), total_timepoints_test)
    expect_true(ncol(result$X_full_design_final) > 0)
    
    # Should have task regressors
    design_colnames <- colnames(result$X_full_design_final)
    expect_true(any(grepl("task_", design_colnames)), "Design matrix should contain task regressors")
  }
  
  # If task betas are extracted, check their structure
  if (!is.null(result$final_task_betas)) {
    expect_equal(ncol(result$final_task_betas), n_voxels_test)
    expect_true(nrow(result$final_task_betas) >= 2)  # At least TaskA and TaskB
  }
})

# Test with minimal but structured data to ensure robustness
test_that("NDX_Process_Subject handles minimal synthetic data gracefully", {
  # Create minimal structured data that should allow basic workflow completion
  set.seed(123)
  n_time_minimal <- 80  # Enough for basic analysis
  n_voxels_minimal <- 10  # More voxels for better analysis
  
  # Generate data with some temporal structure for RPCA to potentially find
  base_signal <- sin(2*pi*seq_len(n_time_minimal) / 20) + 0.5*cos(2*pi*seq_len(n_time_minimal) / 15)
  simple_Y <- matrix(rnorm(n_time_minimal * n_voxels_minimal, sd = 0.5), nrow = n_time_minimal, ncol = n_voxels_minimal)
  
  # Add some structured temporal patterns to a few voxels to give RPCA something to find
  for (v in 1:3) {
    simple_Y[, v] <- simple_Y[, v] + 0.3 * base_signal * runif(1, 0.5, 1.5)
  }
  
  simple_run_idx <- rep(1L, n_time_minimal)
  simple_motion <- matrix(rnorm(n_time_minimal * 3, sd = 0.1), nrow = n_time_minimal, ncol = 3)
  colnames(simple_motion) <- paste0("mot", 1:3)
  
  # Create sufficient events for HRF estimation
  simple_events <- data.frame(
    onsets = c(10, 20, 30, 50, 60, 70) * TR_test,
    durations = rep(4 * TR_test, 6),
    condition = factor(rep(c("TaskA", "TaskB"), 3)),
    blockids = as.integer(rep(1, 6))
  )
  
  simple_options <- list(
    opts_pass0 = list(poly_degree = 1),
    opts_hrf = list(
      hrf_fir_taps = 6, hrf_fir_span_seconds = 12, good_voxel_R2_threshold = -Inf,
      lambda1_grid = c(0.1), lambda2_grid = c(0.1), cv_folds = 2,
      hrf_min_good_voxels = 1, hrf_cluster_method = "none", num_hrf_clusters = 1,
      hrf_min_events_for_fir = 2  # Reasonable threshold for minimal events
    ),
    opts_rpca = list(k_global_target = 2, rpca_lambda_auto = FALSE, rpca_lambda_fixed = 0.01),  # More permissive
    opts_spectral = list(n_sine_candidates = 2, nyquist_guard_factor = 0.8),
    opts_whitening = list(global_ar_on_design = FALSE, max_ar_failures_prop = 0.8),  # Allow more failures
    opts_ridge = list(lambda_ridge = 0.5),
    task_regressor_names_for_extraction = c("task_TaskA", "task_TaskB"),
    max_passes = 1, min_des_gain_convergence = -Inf, min_rho_noise_projection_convergence = -Inf
  )
  
  # Should not error even with minimal data
  expect_no_error({
    result <- NDX_Process_Subject(
      Y_fmri = simple_Y, events = simple_events, motion_params = simple_motion,
      run_idx = simple_run_idx, TR = TR_test, user_options = simple_options, verbose = FALSE
    )
  })
  
  # Basic checks that critical fields are present
  expect_equal(dim(result$Y_residuals_final_unwhitened), dim(simple_Y))
  expect_true(is.numeric(result$pass0_vars) && length(result$pass0_vars) == 1)
  expect_equal(dim(result$pass0_residuals), dim(simple_Y))
  
  # Workflow should complete at least one pass even with minimal data
  expect_true(result$num_passes_completed >= 1)
}) 