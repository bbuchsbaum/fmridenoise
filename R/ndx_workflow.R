#' Process a Single Subject through the Full ND-X Iterative Denoising Workflow
#'
#' This function orchestrates the core modules of the ND-X pipeline.
#' It iteratively performs initial residual generation, FIR HRF estimation,
#' nuisance component identification (RPCA and Spectral), AR(2) pre-whitening,
#' and ridge regression until convergence or max passes. Input types are
#' validated and the function will stop with informative error messages if they
#' do not meet expectations.
#'
#' @param Y_fmri A numeric matrix of fMRI data (total_timepoints x voxels), concatenated across runs if applicable.
#' @param events A data frame describing experimental events. Must contain columns compatible
#'   with `ndx_initial_glm` and `ndx_estimate_initial_hrfs` (e.g., `onsets`, `durations`, `condition`, `blockids`).
#' @param motion_params A numeric matrix of motion parameters (total_timepoints x num_motion_regressors).
#' @param run_idx A numeric vector indicating run membership for each row (timepoint)
#'   in `Y_fmri`, `motion_params`, etc.
#' @param TR Numeric, repetition time in seconds.
#' @param spike_TR_mask Optional. A logical vector of length `nrow(Y_fmri)` where TRUE
#'   indicates a TR to be excluded from HRF estimation and potentially other steps.
#'   If NULL, all TRs are considered valid initially.
#' @param user_options A list containing various sub-lists of user-configurable options for each module:
#'   - `opts_pass0`: List of options for `ndx_initial_glm` (e.g., `poly_degree`). This is also used by `ndx_build_design_matrix` for `poly_degree_val`.
#'   - `opts_hrf`: List of options for `ndx_estimate_initial_hrfs` (e.g., `hrf_fir_taps`, `good_voxel_R2_threshold`, `lambda1_grid`, `lambda2_grid`, `cv_folds`).
#'   - `opts_rpca`: List of options for `ndx_rpca_temporal_components_multirun` (e.g., `k_per_run_target`, `rpca_lambda_auto`).
#'   - `opts_spectral`: List of options for `ndx_spectral_sines` (e.g., `n_sine_candidates`, `nyquist_guard_factor`, `k_tapers`, `nw`).
#'   - `opts_whitening`: List of options for `ndx_ar2_whitening` (e.g., `order`, `global_ar_on_design`, `max_ar_failures_prop`).
#'   - `opts_ridge`: List of options including `lambda_ridge` for `ndx_solve_ridge`.
#'   - `opts_annihilation`: List of options for Annihilation Mode (e.g., `annihilation_enable_mode`, `annihilation_gdlite_k_max`).
#'   - `task_regressor_names_for_extraction` (character vector): Names of task regressors to extract betas for.
#'   - `max_passes` (integer): Maximum number of iterations for the refinement loop. Default: 3.
#'   - `min_des_gain_convergence` (numeric): Minimum DES gain to continue iteration. Default: 0.005.
#'   - `min_rho_noise_projection_convergence` (numeric): Minimum rho noise projection to continue iteration. Default: 0.01. 
#'
#' See `?ndx_user_options` for a complete reference of all available options.
#' @param verbose Logical, if TRUE, print progress messages. Default: TRUE.
#'
#' @return A list containing key outputs from the workflow, such as:
#'   - `final_task_betas`: Extracted task betas after the final pass.
#'   - `diagnostics_per_pass`: A list of diagnostic metrics for each pass.
#'   - `beta_history_per_pass`: A list of beta estimates from each pass.
#'   - (Other final pass outputs like residuals, AR coeffs, nuisance components, etc.)
#'   - `spike_TR_mask`: Logical vector of TRs flagged as spikes from RPCA `S`.
#'   - `pass0_residuals`: Residuals from the initial GLM (Pass 0) for diagnostics.
#'   - `S_matrix_rpca_final`: Final RPCA sparse matrix for spike visualization.
#'   - `ljung_box_p`: Ljung-Box p-value of whitened residuals.
#'   - `num_hrf_clusters`: Effective HRF cluster count from HRF estimation.
#'   - `annihilation_var_ratio` and `annihilation_verdict` when Annihilation mode is used.
#' @importFrom fmrireg event_model design_matrix sampling_frame
#' @importFrom tibble is_tibble
#' @export NDX_Process_Subject
NDX_Process_Subject <- function(Y_fmri,
                              events,
                              motion_params,
                              run_idx,
                              TR,
                              spike_TR_mask = NULL,
                              user_options = list(),
                              verbose = TRUE) {

  if (verbose) message("Starting ND-X Processing for Subject...")

  # Initialize workflow state
  workflow_state <- .ndx_initialize_workflow_state(
    Y_fmri, events, motion_params, run_idx, TR, 
    spike_TR_mask, user_options, verbose
  )
  
  # Run iterative refinement loop
  workflow_state <- .ndx_run_iterative_passes(workflow_state)
  
  # Finalize and return results
  final_results <- .ndx_finalize_workflow_results(workflow_state)
  
  if (verbose) message("ND-X Subject Processing finished.")
  return(final_results)
}

#' Initialize workflow state for NDX processing
#' @keywords internal
.ndx_initialize_workflow_state <- function(Y_fmri, events, motion_params, run_idx, TR,
                                          spike_TR_mask, user_options, verbose) {
  
  ndx_validate_process_subject_inputs(Y_fmri, events, motion_params, run_idx, TR)
  opts <- ndx_prepare_workflow_options(user_options)
  
  # Initialize spike mask
  current_overall_spike_TR_mask <- if (is.null(spike_TR_mask)) {
    rep(FALSE, nrow(Y_fmri))
  } else {
    as.logical(spike_TR_mask)
  }
  
  # Setup annihilation mode if enabled
  annihilation_info <- ndx_run_annihilation_setup(
    Y_fmri, events, motion_params, run_idx, TR,
    opts$opts_annihilation, opts$verbose
  )
  
  # Return initialized state
  list(
    # Input data
    Y_fmri = Y_fmri,
    events = events,
    motion_params = motion_params,
    run_idx = run_idx,
    TR = TR,
    
    # Options
    opts = opts,
    verbose = verbose,
    
    # State variables
    diagnostics_per_pass = list(),
    beta_history_per_pass = list(),
    Y_residuals_current = NULL,
    VAR_BASELINE_FOR_DES = NULL,
    current_overall_spike_TR_mask = current_overall_spike_TR_mask,
    prev_U_NDX_Nuisance = NULL,
    pass0_residuals = NULL,
    
    # Annihilation mode
    U_GD_Lite_fixed_PCs = annihilation_info$selected_pcs,
    gdlite_initial_results = annihilation_info$gdlite_results
  )
}

#' Run iterative refinement passes
#' @keywords internal
.ndx_run_iterative_passes <- function(workflow_state) {
  
  opts <- workflow_state$opts
  max_passes <- opts$max_passes
  verbose <- workflow_state$verbose
  
  for (pass_num in 1:max_passes) {
    if (verbose) message(sprintf("\n--- Starting ND-X Pass %d/%d ---", pass_num, max_passes))
    
    # Run single pass with error handling
    tryCatch({
      workflow_state <- .ndx_run_single_pass(workflow_state, pass_num)
    }, error = function(e) {
      if (verbose) message(sprintf("Pass %d failed with error: %s", pass_num, e$message))
      # Set a minimal current_pass_results if it doesn't exist
      if (is.null(workflow_state$current_pass_results)) {
        workflow_state$current_pass_results <<- list()
      }
      # Break out of the loop by setting Y_residuals_current to NULL
      workflow_state$Y_residuals_current <<- NULL
    })
    
    # Check for early termination
    if (is.null(workflow_state$Y_residuals_current)) {
      if (verbose) message(sprintf("Pass %d: Y_residuals_current became NULL. Stopping iteration.", pass_num))
      break
    }
    
    # Check convergence (only after pass 1)
    if (pass_num > 1) {
      converged <- .ndx_check_convergence(workflow_state, pass_num)
      if (converged) {
        if (verbose) message("Convergence criteria met. Stopping iterations.")
        break
      }
    }
  }
  
  workflow_state$num_passes_completed <- pass_num
  return(workflow_state)
}

#' Run a single pass of the NDX workflow
#' @keywords internal
.ndx_run_single_pass <- function(workflow_state, pass_num) {
  
  current_pass_results <- list()
  
  # Step 1: Initial GLM or use previous residuals
  workflow_state <- .ndx_step_initial_glm(workflow_state, pass_num, current_pass_results)
  
  # Step 2: HRF Estimation
  current_pass_results <- .ndx_step_hrf_estimation(workflow_state, pass_num, current_pass_results)
  
  # Step 3: RPCA Nuisance Components
  result <- .ndx_step_rpca_components(workflow_state, pass_num, current_pass_results)
  workflow_state <- result$workflow_state
  current_pass_results <- result$current_pass_results
  
  # Step 4: Spectral Nuisance Components
  current_pass_results <- .ndx_step_spectral_components(workflow_state, pass_num, current_pass_results)
  
  # Step 5: Construct Design Matrix
  current_pass_results <- .ndx_step_design_matrix(workflow_state, pass_num, current_pass_results)
  
  # Step 6: AR(2) Pre-whitening
  current_pass_results <- .ndx_step_whitening(workflow_state, pass_num, current_pass_results)
  
  # Step 7: Ridge Regression
  result <- .ndx_step_ridge_regression(workflow_state, pass_num, current_pass_results)
  workflow_state <- result$workflow_state
  current_pass_results <- result$current_pass_results
  
  # Step 8: Calculate Diagnostics
  result <- .ndx_step_diagnostics(workflow_state, pass_num, current_pass_results)
  workflow_state <- result$workflow_state
  current_pass_results <- result$current_pass_results
  
  # Store current pass results in workflow state for final output
  workflow_state$current_pass_results <- current_pass_results
  
  return(workflow_state)
}







#' Check convergence criteria
#' @keywords internal
.ndx_check_convergence <- function(workflow_state, pass_num) {
  
  verbose <- workflow_state$verbose
  min_des_gain_convergence <- workflow_state$opts$min_des_gain_convergence
  min_rho_noise_projection_convergence <- workflow_state$opts$min_rho_noise_projection_convergence
  
  converged_des <- FALSE
  converged_rho <- FALSE
  
  # DES gain check
  prev_DES <- workflow_state$diagnostics_per_pass[[pass_num - 1]]$DES
  current_DES <- workflow_state$diagnostics_per_pass[[pass_num]]$DES
  
  if (!is.na(prev_DES) && !is.na(current_DES)) {
    des_gain <- current_DES - prev_DES
    if (verbose) message(sprintf("  Pass %d DES gain: %.4f", pass_num, des_gain))
    if (des_gain < min_des_gain_convergence) {
      if (verbose) message(sprintf("  Convergence potentially met based on DES gain (%.4f < %.4f).", des_gain, min_des_gain_convergence))
      converged_des <- TRUE
    }
  } else {
    if (verbose) message("  Cannot check DES convergence due to NA DES values.")
  }
  
  # Rho noise projection check
  current_rho <- workflow_state$diagnostics_per_pass[[pass_num]]$rho_noise_projection
  if (!is.na(current_rho)) {
    if (verbose) message(sprintf("  Current Rho Noise Projection: %.4f", current_rho))
    if (current_rho < min_rho_noise_projection_convergence) {
      if (verbose) message(sprintf("  Convergence potentially met based on Rho Noise Projection (%.4f < %.4f).", current_rho, min_rho_noise_projection_convergence))
      converged_rho <- TRUE
    }
  } else {
    if (verbose) message("  Cannot check Rho convergence due to NA Rho value.")
  }
  
  return(converged_des || converged_rho)
}

#' Finalize workflow results
#' @keywords internal
.ndx_finalize_workflow_results <- function(workflow_state) {
  
  # Get results from last completed pass
  last_pass_num <- workflow_state$num_passes_completed %||% 0
  current_pass_results <- workflow_state$current_pass_results
  
  # Handle case where workflow failed early
  if (is.null(current_pass_results)) {
    current_pass_results <- list()
  }
  
  # Ensure critical fields are always available with proper fallbacks
  
  # Y_residuals_final_unwhitened: Use current_pass_results first, then workflow_state fallbacks
  final_Y_residuals <- current_pass_results$Y_residuals_final_unwhitened
  if (is.null(final_Y_residuals)) {
    final_Y_residuals <- workflow_state$pass0_residuals
    if (is.null(final_Y_residuals)) {
      final_Y_residuals <- workflow_state$Y_residuals_current
      if (is.null(final_Y_residuals)) {
        # Last resort: use original Y_fmri as residuals (no denoising)
        final_Y_residuals <- workflow_state$Y_fmri
      }
    }
  }
  
  # pass0_residuals: Should always be available from initial GLM
  final_pass0_residuals <- workflow_state$pass0_residuals
  if (is.null(final_pass0_residuals)) {
    final_pass0_residuals <- workflow_state$Y_residuals_current
    if (is.null(final_pass0_residuals)) {
      # Last resort: use original Y_fmri
      final_pass0_residuals <- workflow_state$Y_fmri
    }
  }
  
  # pass0_vars: Should always be available from initial GLM
  final_pass0_vars <- workflow_state$VAR_BASELINE_FOR_DES
  if (is.null(final_pass0_vars)) {
    # Calculate a fallback variance if needed
    if (!is.null(final_pass0_residuals)) {
      final_pass0_vars <- mean(apply(final_pass0_residuals, 2, var, na.rm = TRUE), na.rm = TRUE)
    } else {
      final_pass0_vars <- 1.0  # Default fallback
    }
  }
  
  # Build final results structure with safe defaults
  final_results <- list(
    final_task_betas = current_pass_results$final_task_betas_pass,
    Y_residuals_final_unwhitened = final_Y_residuals,
    ar_coeffs_voxelwise = current_pass_results$ar_coeffs_voxelwise,
    rpca_components = current_pass_results$rpca_components,
    spectral_sines = current_pass_results$spectral_sines,
    estimated_hrfs = current_pass_results$estimated_hrfs,
    S_matrix_rpca_final = current_pass_results$S_matrix_rpca,
    pass0_vars = final_pass0_vars,
    pass0_residuals = final_pass0_residuals,
    na_mask_whitening = current_pass_results$na_mask_whitening,
    spike_TR_mask = workflow_state$current_overall_spike_TR_mask %||% rep(FALSE, nrow(workflow_state$Y_fmri)),
    X_full_design_final = current_pass_results$X_full_design,
    diagnostics_per_pass = workflow_state$diagnostics_per_pass %||% list(),
    beta_history_per_pass = workflow_state$beta_history_per_pass %||% list(),
    num_passes_completed = last_pass_num
  )
  
  # Ensure diagnostics_per_pass has at least one entry if any pass completed
  if (last_pass_num > 0 && length(final_results$diagnostics_per_pass) == 0) {
    # Create a minimal diagnostic entry
    final_results$diagnostics_per_pass <- list(list(
      DES = NA,
      rho_noise_projection = NA,
      ljung_box_p = NA
    ))
  }
  
  # Ensure beta_history_per_pass has at least one entry if any pass completed
  if (last_pass_num > 0 && length(final_results$beta_history_per_pass) == 0) {
    final_results$beta_history_per_pass <- list(NULL)
  }
  
  # Add final diagnostics if available
  if (last_pass_num > 0 && length(workflow_state$diagnostics_per_pass) >= last_pass_num) {
    last_diag <- workflow_state$diagnostics_per_pass[[last_pass_num]]
    final_results$ljung_box_p <- last_diag$ljung_box_p
    final_results$num_hrf_clusters <- last_diag$num_hrf_clusters
  } else {
    # Provide defaults when no diagnostics are available
    final_results$ljung_box_p <- NA
    final_results$num_hrf_clusters <- 1L
  }
  
  # Add annihilation mode results
  final_results <- .ndx_add_annihilation_results(final_results, workflow_state)
  
  return(final_results)
}

#' Add annihilation mode specific results
#' @keywords internal
.ndx_add_annihilation_results <- function(final_results, workflow_state) {
  
  opts_annihilation <- workflow_state$opts$opts_annihilation
  
  if (opts_annihilation$annihilation_enable_mode) {
    final_results$annihilation_mode_active <- TRUE
    final_results$gdlite_pcs <- workflow_state$U_GD_Lite_fixed_PCs
    
    # Add orthogonalized components
    if (!is.null(workflow_state$prev_U_NDX_Nuisance)) {
      final_results$rpca_orthogonalized <- workflow_state$prev_U_NDX_Nuisance$rpca_unique
      final_results$spectral_orthogonalized <- workflow_state$prev_U_NDX_Nuisance$spectral_unique
    }
    
    # Include GLMdenoise diagnostics
    if (!is.null(workflow_state$gdlite_initial_results)) {
      final_results$gdlite_diagnostics <- list(
        noise_pool_mask = workflow_state$gdlite_initial_results$noise_pool_mask,
        good_voxels_mask = workflow_state$gdlite_initial_results$good_voxels_mask,
        tsnr = workflow_state$gdlite_initial_results$tsnr,
        cv_r2 = workflow_state$gdlite_initial_results$cv_r2,
        optimal_k = workflow_state$gdlite_initial_results$optimal_k,
        r2_vals_by_k = workflow_state$gdlite_initial_results$r2_vals_by_k
      )
    }
    
    verdict_stats <- ndx_annihilation_verdict_stats(final_results)
    final_results$annihilation_var_ratio <- verdict_stats$var_ratio
    final_results$annihilation_verdict <- verdict_stats$verdict
  } else {
    final_results$annihilation_mode_active <- FALSE
  }
  
  return(final_results)
}

# Helper for default options, similar to base R's %||%
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}