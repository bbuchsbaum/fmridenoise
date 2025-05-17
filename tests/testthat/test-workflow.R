context("ndx_run_sprint1 - Workflow Integration Tests")

# Mock data and parameters needed for ndx_run_sprint1
# This will be more involved as it requires setting up inputs for several modules.

# Basic setup parameters
TR_test <- 2.0
n_time_per_run_test <- 50 # Shorter for faster tests
n_runs_test <- 1
total_timepoints_test <- n_time_per_run_test * n_runs_test
n_voxels_test <- 5 # Small number of voxels

# Y_fmri
Y_fmri_test <- matrix(rnorm(total_timepoints_test * n_voxels_test), 
                      nrow = total_timepoints_test, 
                      ncol = n_voxels_test)

# run_idx
run_idx_test <- rep(1:n_runs_test, each = n_time_per_run_test)

# motion_params
motion_params_test <- matrix(rnorm(total_timepoints_test * 3), 
                             nrow = total_timepoints_test, ncol = 3)
colnames(motion_params_test) <- paste0("mot", 1:3)

# events data frame
events_test <- data.frame(
  onsets = as.numeric(c(10, 30) * TR_test), # Onsets in seconds
  durations = as.numeric(c(5, 5) * TR_test),   # Durations in seconds
  condition = factor(c("TaskA", "TaskB")),
  blockids = as.integer(rep(1, 2)) # All in the first run
)
if (n_runs_test > 1) {
  events_run2_test <- data.frame(
    onsets = as.numeric(c(10, 30) * TR_test),
    durations = as.numeric(c(5, 5) * TR_test),
    condition = factor(c("TaskA", "TaskB")),
    blockids = as.integer(rep(2, 2))
  )
  events_test <- rbind(events_test, events_run2_test)
}


# User options (minimal for now, can be expanded)
user_options_test <- list(
  opts_pass0 = list(
    poly_degree = 1 # This is used by .construct_final_design_matrix
  ),
  opts_hrf = list(
    hrf_fir_taps = 6,
    hrf_fir_span_seconds = 12, # TR_test * hrf_fir_taps = 2*6=12. Using explicit value.
    good_voxel_R2_threshold = -Inf, # Use all voxels for HRF for simplicity
    lambda1_grid = c(0.1), # Minimal grid for speed
    lambda2_grid = c(0.1),
    cv_folds = 2, # Minimal folds
    hrf_min_good_voxels = 1, # Allow even with 1 voxel for minimal test
    hrf_cluster_method = "none",
    num_hrf_clusters = 1 
  ),
  opts_rpca = list(
    k_global_target = 2, 
    rpca_lambda_auto = FALSE,
    rpca_lambda_fixed = 0.1
  ),
  opts_spectral = list(
    n_sine_candidates = 2, # Small number
    nyquist_guard_factor = 0.1
  ),
  opts_whitening = list(
    global_ar_on_design = FALSE,
    max_ar_failures_prop = 0.5
  ),
  opts_ridge = list(
    lambda_ridge = 0.5
  ),
  task_regressor_names_for_extraction = c("task_TaskA", "task_TaskB"),
  # Workflow control options directly in user_options_test
  max_passes = 2, 
  min_des_gain_convergence = -Inf, 
  min_rho_noise_projection_convergence = -Inf 
)

test_that("NDX_Process_Subject runs with minimal valid inputs and returns correct structure", {
  # skip_on_cran() # Potentially long-running test - REMOVED FOR NOW
  
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
    "final_task_betas", "Y_residuals_final_unwhitened", "ar_coeffs_voxelwise",
    "rpca_components", "spectral_sines", "estimated_hrfs", "S_matrix_rpca_final",
    "pass0_vars", "pass0_residuals", "na_mask_whitening", "spike_TR_mask",
    "X_full_design_final", "diagnostics_per_pass", "beta_history_per_pass",
    "num_passes_completed"
  )
  expect_named(workflow_output, expected = expected_names, ignore.order = TRUE)
  
  # Basic dimension checks for key matrix outputs
  expect_equal(dim(workflow_output$Y_residuals_final_unwhitened), dim(Y_fmri_test))
  expect_true(is.numeric(workflow_output$pass0_vars) && length(workflow_output$pass0_vars) == 1)
  expect_equal(dim(workflow_output$pass0_residuals), dim(Y_fmri_test))
  
  if (!is.null(workflow_output$X_full_design_final)) {
      expect_equal(nrow(workflow_output$X_full_design_final), total_timepoints_test)
  }

  expect_length(workflow_output$spike_TR_mask, total_timepoints_test)
  expect_type(workflow_output$spike_TR_mask, "logical")
  
  # Check diagnostics_per_pass structure
  expect_true(is.list(workflow_output$diagnostics_per_pass))
  expect_length(workflow_output$diagnostics_per_pass, workflow_output$num_passes_completed)
  if (workflow_output$num_passes_completed > 0) {
    expect_true(all(sapply(workflow_output$diagnostics_per_pass, function(p) "DES" %in% names(p))))
    if (workflow_output$num_passes_completed > 1) { # Rho is calculated from pass 1 onwards effectively for convergence
         # Rho might be NA if no nuisance components or residuals in a pass
         # expect_true(all(sapply(workflow_output$diagnostics_per_pass, function(p) "rho_noise_projection" %in% names(p))))
    }
  }
  expect_true(is.list(workflow_output$beta_history_per_pass))
  expect_length(workflow_output$beta_history_per_pass, workflow_output$num_passes_completed)
  
  # Check dimensions of whitened data if X_full_design_final was created
  # Note: Y_whitened and X_whitened are intermediate from the last pass, not primary outputs now
  # We can check their presence in the last element of a more detailed per-pass result if needed later.
  
  # Check ar_coeffs if present
  if (!is.null(workflow_output$ar_coeffs_voxelwise)) {
      expect_true(is.matrix(workflow_output$ar_coeffs_voxelwise))
      expect_equal(ncol(workflow_output$ar_coeffs_voxelwise), n_voxels_test)
      # Expect 2 rows for AR(2) coefficients, or 3 if intercept is included by ar.yw
      # ndx_ar2_whitening stores only the 2 AR coefficients. 
      expect_equal(nrow(workflow_output$ar_coeffs_voxelwise), 2)
  }

  if (!is.null(workflow_output$S_matrix_rpca_final)) {
      expect_equal(dim(workflow_output$S_matrix_rpca_final), dim(Y_fmri_test))
  }

  # Add more specific checks as needed, e.g., for dimensions of rpca/spectral components,
  # and type/structure of estimated_hrfs and final_task_betas.
  # For now, primarily checking it runs and returns expected named elements.
})

# Add more tests: e.g. with spike_TR_mask, different option combinations, edge cases. 