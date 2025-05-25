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

#' Step 1: Initial GLM or use previous residuals
#' @keywords internal
.ndx_step_initial_glm <- function(workflow_state, pass_num, current_pass_results) {
  
  verbose <- workflow_state$verbose
  
  if (pass_num == 1) {
    if (verbose) message(sprintf("Pass %d: Running Initial GLM for residual generation...", pass_num))
    
    initial_glm_output <- ndx_initial_glm(
      Y_fmri = workflow_state$Y_fmri,
      events = workflow_state$events,
      motion_params = workflow_state$motion_params,
      run_idx = workflow_state$run_idx,
      TR = workflow_state$TR,
      poly_degree = workflow_state$opts$opts_pass0$poly_degree
    )
    
    workflow_state$Y_residuals_current <- initial_glm_output$Y_residuals_current
    workflow_state$VAR_BASELINE_FOR_DES <- initial_glm_output$VAR_BASELINE_FOR_DES
    workflow_state$pass0_residuals <- initial_glm_output$Y_residuals_current
    
    current_pass_results$Y_residuals_from_glm <- workflow_state$Y_residuals_current
    current_pass_results$pass0_vars <- workflow_state$VAR_BASELINE_FOR_DES
    
  } else {
    if (verbose) message(sprintf("Pass %d: Using residuals from Pass %d.", pass_num, pass_num - 1))
    
    if (is.null(workflow_state$Y_residuals_current)) {
      warning(sprintf("Pass %d: Y_residuals_current is NULL. Cannot proceed. Stopping iterations.", pass_num))
      return(workflow_state)
    }
  }
  
  return(workflow_state)
}

#' Step 2: HRF Estimation
#' @keywords internal
.ndx_step_hrf_estimation <- function(workflow_state, pass_num, current_pass_results) {
  
  verbose <- workflow_state$verbose
  
  if (verbose) message(sprintf("Pass %d: Estimating FIR HRFs...", pass_num))
  
  estimated_hrfs_raw <- ndx_estimate_initial_hrfs(
    Y_fmri = workflow_state$Y_fmri,
    pass0_residuals = workflow_state$Y_residuals_current,
    events = workflow_state$events,
    run_idx = workflow_state$run_idx,
    TR = workflow_state$TR,
    spike_TR_mask = workflow_state$current_overall_spike_TR_mask,
    user_options = workflow_state$opts$opts_hrf
  )
  
  current_pass_results$estimated_hrfs_per_cluster <- estimated_hrfs_raw
  current_pass_results$num_hrf_clusters <- attr(estimated_hrfs_raw, "num_effective_clusters") %||% 1L
  
  # Process HRFs for design matrix
  current_pass_results$estimated_hrfs <- .ndx_process_hrfs_for_design(
    estimated_hrfs_raw, workflow_state$events, pass_num, verbose
  )
  
  return(current_pass_results)
}

#' Process HRFs for design matrix construction
#' @keywords internal
.ndx_process_hrfs_for_design <- function(estimated_hrfs_raw, events, pass_num, verbose) {
  
  if (!is.null(estimated_hrfs_raw) && tibble::is_tibble(estimated_hrfs_raw) && nrow(estimated_hrfs_raw) > 0) {
    if (verbose) message(sprintf("  Pass %d: Received %d raw HRF estimates (per cluster/condition).", pass_num, nrow(estimated_hrfs_raw)))
    
    if (any(colnames(estimated_hrfs_raw) == "cluster_id")) {
      # Use cluster 1's HRF for each condition
      estimated_hrfs_for_design <- estimated_hrfs_raw[
        estimated_hrfs_raw$cluster_id == 1L, 
        c("condition", "hrf_estimate", "taps")
      ]
      estimated_hrfs_for_design <- unique(tibble::as_tibble(estimated_hrfs_for_design))
    } else {
      # No clustering, use all HRFs
      estimated_hrfs_for_design <- tibble::as_tibble(
        estimated_hrfs_raw[, c("condition", "hrf_estimate", "taps")]
      )
    }
    
    # Check for missing HRFs
    original_conditions <- unique(as.character(events$condition))
    missing_hrfs <- original_conditions[!original_conditions %in% estimated_hrfs_for_design$condition]
    if (length(missing_hrfs) > 0 && verbose) {
      message(sprintf("    Pass %d: Not all conditions have a selected HRF for design matrix. Missing: %s", 
                     pass_num, paste(missing_hrfs, collapse=", ")))
    }
    
    return(estimated_hrfs_for_design)
  } else {
    if (verbose) message(sprintf("Pass %d: HRF estimation returned NULL or empty. No task regressors will be built for design matrix.", pass_num))
    return(NULL)
  }
}

#' Step 3: RPCA Nuisance Components
#' @keywords internal
.ndx_step_rpca_components <- function(workflow_state, pass_num, current_pass_results) {
  
  verbose <- workflow_state$verbose
  opts_rpca <- workflow_state$opts$opts_rpca
  
  if (verbose) message(sprintf("Pass %d: Identifying RPCA nuisance components...", pass_num))
  
  # Prepare adaptive RPCA options
  current_opts_rpca <- .ndx_prepare_adaptive_rpca_options(workflow_state, pass_num, opts_rpca)
  
  # Determine k_rpca_global
  k_rpca_global <- .ndx_determine_rpca_k(workflow_state, pass_num, current_opts_rpca)
  current_pass_results$k_rpca_global <- k_rpca_global
  
  # Run RPCA
  rpca_out <- ndx_rpca_temporal_components_multirun(
    Y_residuals_cat = workflow_state$Y_residuals_current,
    run_idx = workflow_state$run_idx,
    k_global_target = k_rpca_global,
    user_options = current_opts_rpca
  )
  
  # Process RPCA results
  workflow_state <- .ndx_process_rpca_results(workflow_state, rpca_out, current_pass_results, verbose)
  
  return(list(workflow_state = workflow_state, current_pass_results = current_pass_results))
}

#' Prepare adaptive RPCA options
#' @keywords internal
.ndx_prepare_adaptive_rpca_options <- function(workflow_state, pass_num, opts_rpca) {
  
  current_opts_rpca <- opts_rpca
  verbose <- workflow_state$verbose
  
  if (pass_num > 1 && !is.null(workflow_state$diagnostics_per_pass[[pass_num - 1]]$V_global_singular_values_from_rpca)) {
    current_opts_rpca$adaptive_k_singular_values <- workflow_state$diagnostics_per_pass[[pass_num - 1]]$V_global_singular_values_from_rpca
    if (verbose) {
      message(sprintf("    Pass %d: Using %d singular values from previous pass for RPCA rank adaptation.", 
                     pass_num, length(current_opts_rpca$adaptive_k_singular_values)))
    }
  } else if (pass_num > 1 && verbose) {
    message(sprintf("    Pass %d: No singular values from previous pass available for RPCA rank adaptation.", pass_num))
  }
  
  return(current_opts_rpca)
}

#' Determine RPCA k value
#' @keywords internal
.ndx_determine_rpca_k <- function(workflow_state, pass_num, current_opts_rpca) {
  
  verbose <- workflow_state$verbose
  
  if (pass_num > 1 && !is.null(current_opts_rpca$adaptive_k_singular_values)) {
    sv_for_adapt <- current_opts_rpca$adaptive_k_singular_values
    if (is.numeric(sv_for_adapt) && length(sv_for_adapt) > 0) {
      drop_ratio <- current_opts_rpca$k_elbow_drop_ratio %||% 0.02
      k_min <- current_opts_rpca$k_rpca_min %||% 20L
      k_max <- current_opts_rpca$k_rpca_max %||% 50L
      
      k_rpca_global <- Auto_Adapt_RPCA_Rank(
        sv_for_adapt, 
        drop_ratio = drop_ratio,
        k_min = k_min,
        k_max = k_max
      )
      if (verbose) message(sprintf("    Pass %d: Adaptive k_rpca_global set to %d", pass_num, k_rpca_global))
    } else {
      k_rpca_global <- current_opts_rpca$k_global_target %||% 5
      if (verbose) message(sprintf("    Pass %d: Invalid/NULL singular values for adaptation, using k_rpca_global = %d from opts/default", pass_num, k_rpca_global))
    }
  } else {
    k_rpca_global <- current_opts_rpca$k_global_target %||% 5
    if (verbose && pass_num == 1) {
      message(sprintf("    Pass %d: Using initial k_rpca_global = %d (from opts/default)", pass_num, k_rpca_global))
    } else if (verbose && pass_num > 1) {
      message(sprintf("    Pass %d: No adaptive singular values provided, using k_rpca_global = %d (from opts/default)", pass_num, k_rpca_global))
    }
  }
  
  return(k_rpca_global)
}

#' Process RPCA results
#' @keywords internal
.ndx_process_rpca_results <- function(workflow_state, rpca_out, current_pass_results, verbose) {
  
  if (!is.null(rpca_out) && is.list(rpca_out)) {
    rpca_components <- rpca_out$C_components
    
    # Update spike mask
    if (!is.null(rpca_out$spike_TR_mask) && is.logical(rpca_out$spike_TR_mask) && 
        length(rpca_out$spike_TR_mask) == length(workflow_state$current_overall_spike_TR_mask)) {
      workflow_state$current_overall_spike_TR_mask <- workflow_state$current_overall_spike_TR_mask | rpca_out$spike_TR_mask
    } else if (!is.null(rpca_out$spike_TR_mask) && verbose) {
      message("    Warning: spike_TR_mask from rpca_out is not a valid logical vector of correct length. Global spike mask not updated by this pass's RPCA.")
    }
    
    current_pass_results$S_matrix_rpca <- rpca_out$S_matrix_cat
    current_pass_results$V_global_singular_values_from_rpca <- rpca_out$V_global_singular_values
    
    # Process precision weights
    if (!is.null(rpca_out$S_matrix_cat) && is.matrix(rpca_out$S_matrix_cat) && is.numeric(rpca_out$S_matrix_cat)) {
      current_pass_results$precision_weights <- ndx_precision_weights_from_S(rpca_out$S_matrix_cat)
    } else {
      if (verbose) message("    RPCA did not return a valid numeric S matrix; skipping precision weighting for this pass.")
      current_pass_results$precision_weights <- NULL
    }
  } else {
    if (verbose && !is.null(rpca_out)) {
      message("    Warning: rpca_out from ndx_rpca_temporal_components_multirun was not NULL but not a list as expected.")
    } else if (verbose && is.null(rpca_out)) {
      message("    ndx_rpca_temporal_components_multirun returned NULL (RPCA likely failed or yielded no components).")
    }
    rpca_components <- NULL
    current_pass_results$S_matrix_rpca <- NULL
    current_pass_results$V_global_singular_values_from_rpca <- NULL
    current_pass_results$precision_weights <- NULL
  }
  
  current_pass_results$rpca_components <- rpca_components
  current_pass_results$spike_TR_mask <- workflow_state$current_overall_spike_TR_mask
  
  return(workflow_state)
}

#' Step 4: Spectral Nuisance Components
#' @keywords internal
.ndx_step_spectral_components <- function(workflow_state, pass_num, current_pass_results) {
  
  verbose <- workflow_state$verbose
  opts_spectral <- workflow_state$opts$opts_spectral
  
  if (verbose) message(sprintf("Pass %d: Identifying spectral nuisance components...", pass_num))
  
  if (!is.null(workflow_state$Y_residuals_current) && ncol(workflow_state$Y_residuals_current) > 0) {
    mean_residual_for_spectrum <- rowMeans(workflow_state$Y_residuals_current, na.rm = TRUE)
    spectral_sines <- ndx_spectral_sines(
      mean_residual_for_spectrum = mean_residual_for_spectrum,
      TR = workflow_state$TR,
      n_sine_candidates = opts_spectral$n_sine_candidates %||% 10,
      nyquist_guard_factor = opts_spectral$nyquist_guard_factor %||% 0.9,
      k_tapers = opts_spectral$k_tapers %||% 5,
      nw = opts_spectral$nw %||% 3
    )
    current_pass_results$spectral_sines <- spectral_sines
    current_pass_results$num_spectral_sines <- if (!is.null(spectral_sines)) ncol(spectral_sines)/2 else 0
  } else {
    if (verbose) message(sprintf("  Pass %d: Skipping spectral analysis due to no/empty residuals.", pass_num))
    current_pass_results$spectral_sines <- NULL
  }
  
  return(current_pass_results)
}

#' Step 5: Construct Design Matrix
#' @keywords internal
.ndx_step_design_matrix <- function(workflow_state, pass_num, current_pass_results) {
  
  verbose <- workflow_state$verbose
  opts_annihilation <- workflow_state$opts$opts_annihilation
  
  if (verbose) message(sprintf("Pass %d: Constructing full design matrix...", pass_num))
  
  if (opts_annihilation$annihilation_enable_mode && !is.null(workflow_state$U_GD_Lite_fixed_PCs) && 
      ncol(workflow_state$U_GD_Lite_fixed_PCs) > 0) {
    
    X_full_design <- .ndx_build_annihilation_design_matrix(workflow_state, current_pass_results, verbose)
  } else {
    X_full_design <- .ndx_build_standard_design_matrix(workflow_state, current_pass_results)
  }
  
  current_pass_results$X_full_design <- X_full_design
  return(current_pass_results)
}

#' Build design matrix for annihilation mode
#' @keywords internal
.ndx_build_annihilation_design_matrix <- function(workflow_state, current_pass_results, verbose) {
  
  # Standard components for base design
  X_full_design <- ndx_build_design_matrix(
    estimated_hrfs = current_pass_results$estimated_hrfs,
    events = workflow_state$events,
    motion_params = workflow_state$motion_params,
    rpca_components = NULL,  # Add separately with different prefixes
    spectral_sines = NULL,   # Add separately with different prefixes
    run_idx = workflow_state$run_idx,
    TR = workflow_state$TR,
    poly_degree = workflow_state$opts$opts_pass0$poly_degree,
    verbose = FALSE
  )
  
  if (!is.null(X_full_design)) {
    # Add GLMdenoise PCs
    if (ncol(workflow_state$U_GD_Lite_fixed_PCs) > 0) {
      gdlite_pcols <- ncol(workflow_state$U_GD_Lite_fixed_PCs)
      colnames_gdlite <- paste0("gdlite_pc_", seq_len(gdlite_pcols))
      X_gdlite <- workflow_state$U_GD_Lite_fixed_PCs
      colnames(X_gdlite) <- colnames_gdlite
      X_full_design <- cbind(X_full_design, X_gdlite)
      if (verbose) message(sprintf("    Added %d GLMdenoise PCs to design matrix.", gdlite_pcols))
    }
    
    # Add orthogonalized components
    X_full_design <- .ndx_add_orthogonalized_components(X_full_design, workflow_state, verbose)
  }
  
  return(X_full_design)
}

#' Add orthogonalized components to design matrix
#' @keywords internal
.ndx_add_orthogonalized_components <- function(X_full_design, workflow_state, verbose) {
  
  U_NDX_Nuisance_Combined_list <- workflow_state$prev_U_NDX_Nuisance
  if (is.null(U_NDX_Nuisance_Combined_list)) U_NDX_Nuisance_Combined_list <- list()
  
  # Add orthogonalized RPCA components
  if (!is.null(U_NDX_Nuisance_Combined_list$rpca_unique) && 
      ncol(U_NDX_Nuisance_Combined_list$rpca_unique) > 0) {
    
    rpca_unique_pcols <- ncol(U_NDX_Nuisance_Combined_list$rpca_unique)
    colnames_rpca_unique <- paste0("rpca_unique_comp_", seq_len(rpca_unique_pcols))
    X_rpca_unique <- U_NDX_Nuisance_Combined_list$rpca_unique
    colnames(X_rpca_unique) <- colnames_rpca_unique
    X_full_design <- cbind(X_full_design, X_rpca_unique)
    if (verbose) message(sprintf("    Added %d orthogonalized RPCA components to design matrix.", rpca_unique_pcols))
  }
  
  # Add orthogonalized spectral components
  if (!is.null(U_NDX_Nuisance_Combined_list$spectral_unique) && 
      ncol(U_NDX_Nuisance_Combined_list$spectral_unique) > 0) {
    
    spectral_unique_pcols <- ncol(U_NDX_Nuisance_Combined_list$spectral_unique)
    colnames_spectral_unique <- paste0("spectral_unique_comp_", seq_len(spectral_unique_pcols))
    X_spectral_unique <- U_NDX_Nuisance_Combined_list$spectral_unique
    colnames(X_spectral_unique) <- colnames_spectral_unique
    X_full_design <- cbind(X_full_design, X_spectral_unique)
    if (verbose) message(sprintf("    Added %d orthogonalized spectral components to design matrix.", spectral_unique_pcols))
  }
  
  return(X_full_design)
}

#' Build standard design matrix
#' @keywords internal
.ndx_build_standard_design_matrix <- function(workflow_state, current_pass_results) {
  
  X_full_design <- ndx_build_design_matrix(
    estimated_hrfs = current_pass_results$estimated_hrfs,
    events = workflow_state$events,
    motion_params = workflow_state$motion_params,
    rpca_components = current_pass_results$rpca_components,
    spectral_sines = current_pass_results$spectral_sines,
    run_idx = workflow_state$run_idx,
    TR = workflow_state$TR,
    poly_degree = workflow_state$opts$opts_pass0$poly_degree,
    verbose = FALSE
  )
  
  return(X_full_design)
}

#' Step 6: AR(2) Pre-whitening
#' @keywords internal
.ndx_step_whitening <- function(workflow_state, pass_num, current_pass_results) {
  
  verbose <- workflow_state$verbose
  opts_whitening <- workflow_state$opts$opts_whitening
  
  if (verbose) message(sprintf("Pass %d: Performing AR(2) pre-whitening...", pass_num))
  
  if (!is.null(current_pass_results$X_full_design)) {
    # Get residuals for AR fit
    temp_glm_for_ar_residuals <- tryCatch({
      stats::lm.fit(current_pass_results$X_full_design, workflow_state$Y_fmri)$residuals
    }, error = function(e) {
      if(verbose) message(sprintf("  Pass %d: Error fitting temporary GLM for AR residuals: %s. Using Y_residuals_current for AR fit.", pass_num, e$message))
      workflow_state$Y_residuals_current
    })
    
    whitening_output <- ndx_ar2_whitening(
      Y_data = workflow_state$Y_fmri,
      X_design_full = current_pass_results$X_full_design,
      Y_residuals_for_AR_fit = temp_glm_for_ar_residuals,
      order = opts_whitening$order %||% 2L,
      global_ar_on_design = opts_whitening$global_ar_on_design %||% TRUE,
      weights = NULL,  # Remove precision weights from AR estimation - they're most effective in ridge regression
      max_ar_failures_prop = opts_whitening$max_ar_failures_prop %||% 0.3
    )
    
    current_pass_results$Y_whitened <- whitening_output$Y_whitened
    current_pass_results$X_whitened <- whitening_output$X_whitened
    current_pass_results$ar_coeffs_voxelwise <- whitening_output$ar_coeffs_voxelwise
    current_pass_results$na_mask_whitening <- whitening_output$na_mask
  } else {
    if (verbose) message(sprintf("  Pass %d: Skipping AR(2) whitening as X_full_design is NULL.", pass_num))
    current_pass_results$Y_whitened <- workflow_state$Y_fmri
    current_pass_results$X_whitened <- current_pass_results$X_full_design
    current_pass_results$ar_coeffs_voxelwise <- NULL
    current_pass_results$na_mask_whitening <- rep(FALSE, nrow(workflow_state$Y_fmri))
  }
  
  return(current_pass_results)
}

#' Step 7: Ridge Regression
#' @keywords internal
.ndx_step_ridge_regression <- function(workflow_state, pass_num, current_pass_results) {
  
  verbose <- workflow_state$verbose
  opts_ridge <- workflow_state$opts$opts_ridge
  task_regressor_names <- workflow_state$opts$task_regressor_names_for_extraction
  
  if (verbose) message(sprintf("Pass %d: Performing Ridge/Anisotropic Regression...", pass_num))
  
  if (!is.null(current_pass_results$X_whitened) && !is.null(current_pass_results$Y_whitened) && 
      ncol(current_pass_results$X_whitened) > 0) {
    
    # Run ridge regression
    ridge_results <- .ndx_run_ridge_regression(current_pass_results, opts_ridge, verbose)
    current_pass_results$ridge_betas_whitened <- ridge_results$ridge_betas_whitened
    current_pass_results$lambda_parallel_noise <- ridge_results$lambda_parallel_noise
    current_pass_results$lambda_perp_signal <- ridge_results$lambda_perp_signal
    
    # Store beta history
    workflow_state$beta_history_per_pass[[pass_num]] <- ridge_results$ridge_betas_whitened
    
    # Extract task betas
    if (!is.null(ridge_results$ridge_betas_whitened) && length(task_regressor_names) > 0) {
      final_task_betas_pass <- ndx_extract_task_betas(
        betas_whitened = ridge_results$ridge_betas_whitened,
        X_whitened_colnames = colnames(current_pass_results$X_whitened), 
        task_regressor_names = task_regressor_names,
        ar_coeffs_global = NULL 
      )
      current_pass_results$final_task_betas_pass <- final_task_betas_pass
    } else {
      current_pass_results$final_task_betas_pass <- NULL
    }
    
    # Calculate residuals for next pass
    workflow_state <- .ndx_calculate_pass_residuals(workflow_state, current_pass_results)
    
  } else {
    if (verbose) message(sprintf("  Pass %d: Skipping Ridge Regression due to missing whitened data/design.", pass_num))
    current_pass_results$ridge_betas_whitened <- NULL
    current_pass_results$final_task_betas_pass <- NULL
    workflow_state$beta_history_per_pass[[pass_num]] <- NULL
    workflow_state$Y_residuals_current <- NULL
    current_pass_results$Y_residuals_final_unwhitened <- workflow_state$Y_residuals_current
  }
  
  return(list(workflow_state = workflow_state, current_pass_results = current_pass_results))
}

#' Run ridge regression with anisotropic or isotropic penalties
#' @keywords internal
.ndx_run_ridge_regression <- function(current_pass_results, opts_ridge, verbose) {
  
  n_regressors <- ncol(current_pass_results$X_whitened)
  precision_weights_for_pass <- current_pass_results$precision_weights
  
  if (opts_ridge$anisotropic_ridge_enable %||% TRUE) {
    ridge_results <- .ndx_run_anisotropic_ridge(current_pass_results, opts_ridge, precision_weights_for_pass, verbose)
  } else {
    ridge_results <- .ndx_run_isotropic_ridge(current_pass_results, opts_ridge, precision_weights_for_pass, n_regressors)
  }
  
  return(ridge_results)
}

#' Run anisotropic ridge regression
#' @keywords internal
.ndx_run_anisotropic_ridge <- function(current_pass_results, opts_ridge, precision_weights_for_pass, verbose) {
  
  regressor_names <- colnames(current_pass_results$X_whitened)
  n_regressors <- ncol(current_pass_results$X_whitened)
  
  if (!is.null(regressor_names)) {
    # Identify regressor types
    regressor_indices <- .ndx_identify_regressor_types(regressor_names)
    
    # Build projection matrices
    proj_mats <- .ndx_build_projection_matrices(regressor_indices, n_regressors, opts_ridge)
    
    # Estimate residual variance
    res_var_est <- ndx_estimate_res_var_whitened(
      current_pass_results$Y_whitened,
      current_pass_results$X_whitened,
      current_pass_results$na_mask_whitening
    )
    
    # Tune lambda parameters
    lambda_results <- .ndx_tune_lambda_parameters(current_pass_results, proj_mats, opts_ridge, res_var_est)
    
    # Solve anisotropic ridge
    ridge_res <- ndx_solve_anisotropic_ridge(
      Y_whitened = current_pass_results$Y_whitened,
      X_whitened = current_pass_results$X_whitened,
      projection_mats = proj_mats,
      lambda_values = lambda_results$lambda_values,
      na_mask = current_pass_results$na_mask_whitening,
      weights = precision_weights_for_pass,
      gcv_lambda = FALSE,
      res_var_scale = res_var_est
    )
    ridge_betas_whitened <- ridge_res$betas

    return(list(
      ridge_betas_whitened = ridge_betas_whitened,
      lambda_parallel_noise = lambda_results$lambda_parallel_tuned,
      lambda_perp_signal = lambda_results$lambda_perp_signal
    ))
  } else {
    if (verbose) message("    Colnames for X_whitened not available, anisotropic ridge skipped.")
    return(list(ridge_betas_whitened = NULL, lambda_parallel_noise = NULL, lambda_perp_signal = NULL))
  }
}

#' Identify regressor types for anisotropic ridge
#' @keywords internal
.ndx_identify_regressor_types <- function(regressor_names) {
  
  is_rpca_col <- grepl("^rpca_comp_", regressor_names)
  is_spectral_col <- grepl("^(sin_f|cos_f|spectral_comp_)", regressor_names)
  is_gdlite_col <- grepl("^gdlite_pc_", regressor_names)
  is_rpca_unique_col <- grepl("^rpca_unique_comp_", regressor_names)
  is_spectral_unique_col <- grepl("^spectral_unique_comp_", regressor_names)
  
  list(
    noise_rpca_col_indices = which(is_rpca_col),
    noise_spectral_col_indices = which(is_spectral_col),
    noise_gdlite_col_indices = which(is_gdlite_col),
    noise_rpca_unique_col_indices = which(is_rpca_unique_col),
    noise_spectral_unique_col_indices = which(is_spectral_unique_col)
  )
}

#' Build projection matrices for anisotropic ridge
#' @keywords internal
.ndx_build_projection_matrices <- function(regressor_indices, n_regressors, opts_ridge) {
  
  I_reg <- diag(1, n_regressors)
  
  # Build U_noise
  U_noise <- NULL
  if (length(regressor_indices$noise_rpca_col_indices) > 0)
    U_noise <- cbind(U_noise, I_reg[, regressor_indices$noise_rpca_col_indices, drop = FALSE])
  if (length(regressor_indices$noise_spectral_col_indices) > 0)
    U_noise <- cbind(U_noise, I_reg[, regressor_indices$noise_spectral_col_indices, drop = FALSE])
  
  # Build U_gd
  U_gd <- NULL
  if (opts_ridge$annihilation_enable_mode && length(regressor_indices$noise_gdlite_col_indices) > 0)
    U_gd <- I_reg[, regressor_indices$noise_gdlite_col_indices, drop = FALSE]
  
  # Build U_unique
  U_unique <- NULL
  if (opts_ridge$annihilation_enable_mode) {
    if (length(regressor_indices$noise_rpca_unique_col_indices) > 0)
      U_unique <- cbind(U_unique, I_reg[, regressor_indices$noise_rpca_unique_col_indices, drop = FALSE])
    if (length(regressor_indices$noise_spectral_unique_col_indices) > 0)
      U_unique <- cbind(U_unique, I_reg[, regressor_indices$noise_spectral_unique_col_indices, drop = FALSE])
  }
  
  ndx_compute_projection_matrices(U_GD = U_gd, U_Unique = U_unique, U_Noise = U_noise, n_regressors = n_regressors)
}

#' Tune lambda parameters for anisotropic ridge
#' @keywords internal
.ndx_tune_lambda_parameters <- function(current_pass_results, proj_mats, opts_ridge, res_var_est) {
  
  lambda_ratio <- (opts_ridge$lambda_perp_signal %||% 0.1) / (opts_ridge$lambda_parallel_noise %||% 10.0)
  lambda_grid <- 10^seq(-2, 2, length.out = 5) * (res_var_est %||% 1)
  
  lambda_parallel_tuned <- ndx_gcv_tune_lambda_parallel(
    current_pass_results$Y_whitened,
    current_pass_results$X_whitened,
    proj_mats$P_Noise,
    lambda_grid = lambda_grid,
    lambda_ratio = lambda_ratio
  )
  
  lambda_values <- list(
    lambda_parallel = lambda_parallel_tuned,
    lambda_perp_signal = lambda_ratio * lambda_parallel_tuned,
    lambda_gd = opts_ridge$lambda_noise_gdlite %||% lambda_parallel_tuned,
    lambda_unique = opts_ridge$lambda_noise_ndx_unique %||% lambda_parallel_tuned
  )
  
  list(
    lambda_values = lambda_values,
    lambda_parallel_tuned = lambda_parallel_tuned,
    lambda_perp_signal = lambda_ratio * lambda_parallel_tuned
  )
}

#' Run isotropic ridge regression
#' @keywords internal
.ndx_run_isotropic_ridge <- function(current_pass_results, opts_ridge, precision_weights_for_pass, n_regressors) {
  
  K_penalty_diag <- rep(opts_ridge$lambda_ridge %||% 1.0, n_regressors)
  lambda_value <- opts_ridge$lambda_ridge %||% 1.0
  
  ridge_res <- ndx_solve_anisotropic_ridge(
    Y_whitened = current_pass_results$Y_whitened,
    X_whitened = current_pass_results$X_whitened,
    K_penalty_diag = K_penalty_diag,
    na_mask = current_pass_results$na_mask_whitening,
    weights = precision_weights_for_pass
  )

  list(
    ridge_betas_whitened = ridge_res$betas,
    lambda_parallel_noise = lambda_value,
    lambda_perp_signal = lambda_value
  )
}

#' Calculate residuals for next pass
#' @keywords internal
.ndx_calculate_pass_residuals <- function(workflow_state, current_pass_results) {
  
  if(!is.null(current_pass_results$X_full_design) && !is.null(current_pass_results$ridge_betas_whitened)) {
    Y_predicted_unwhitened <- current_pass_results$X_full_design %*% current_pass_results$ridge_betas_whitened
    workflow_state$Y_residuals_current <- workflow_state$Y_fmri - Y_predicted_unwhitened
    current_pass_results$Y_residuals_final_unwhitened <- workflow_state$Y_residuals_current
  } else {
    workflow_state$Y_residuals_current <- NULL
    current_pass_results$Y_residuals_final_unwhitened <- workflow_state$Y_residuals_current
  }
  
  return(workflow_state)
}

#' Step 8: Calculate Diagnostics
#' @keywords internal
.ndx_step_diagnostics <- function(workflow_state, pass_num, current_pass_results) {
  
  verbose <- workflow_state$verbose
  
  if (verbose) message(sprintf("Pass %d: Calculating diagnostics...", pass_num))
  
  pass_diagnostics <- list()
  
  # Calculate Ljung-Box test
  pass_diagnostics$ljung_box_p <- .ndx_calculate_ljung_box(current_pass_results)
  
  # Calculate DES
  pass_diagnostics$DES <- .ndx_calculate_des(workflow_state, pass_num, verbose)
  
  # Calculate Rho Noise Projection and handle annihilation mode
  rho_results <- .ndx_calculate_rho_noise_projection(workflow_state, current_pass_results, verbose)
  pass_diagnostics$rho_noise_projection <- rho_results$rho_noise_projection
  annihilation_active_this_pass <- rho_results$annihilation_active_this_pass
  current_pass_results$U_NDX_Nuisance_Combined_list <- rho_results$U_NDX_Nuisance_Combined_list
  
  # Store other diagnostics
  pass_diagnostics <- .ndx_store_additional_diagnostics(pass_diagnostics, current_pass_results, annihilation_active_this_pass)
  
  # Store diagnostics and update workflow state
  workflow_state$diagnostics_per_pass[[pass_num]] <- pass_diagnostics
  workflow_state$prev_U_NDX_Nuisance <- current_pass_results$U_NDX_Nuisance_Combined_list
  
  return(list(workflow_state = workflow_state, current_pass_results = current_pass_results))
}

#' Calculate Ljung-Box test statistic
#' @keywords internal
.ndx_calculate_ljung_box <- function(current_pass_results) {
  
  whitened_resids <- NULL
  if (!is.null(current_pass_results$ridge_betas_whitened) &&
      !is.null(current_pass_results$Y_whitened) &&
      !is.null(current_pass_results$X_whitened)) {
    
    whitened_resids <- current_pass_results$Y_whitened - current_pass_results$X_whitened %*% current_pass_results$ridge_betas_whitened
    if (!is.null(current_pass_results$na_mask_whitening)) {
      whitened_resids <- whitened_resids[!current_pass_results$na_mask_whitening, , drop = FALSE]
    }
  }
  
  ndx_ljung_box_pval(whitened_resids)
}

#' Calculate DES (Denoising Efficacy Score)
#' @keywords internal
.ndx_calculate_des <- function(workflow_state, pass_num, verbose) {
  
  if (!is.null(workflow_state$Y_residuals_current) && !is.null(workflow_state$VAR_BASELINE_FOR_DES)) {
    des_value <- calculate_DES(
      current_residuals_unwhitened = workflow_state$Y_residuals_current, 
      VAR_BASELINE_FOR_DES = workflow_state$VAR_BASELINE_FOR_DES
    )
    if (verbose) message(sprintf("  Pass %d DES: %.4f", pass_num, des_value))
    return(des_value)
  } else {
    if (verbose) message(sprintf("  Pass %d DES: NA (missing residuals or baseline variance)", pass_num))
    return(NA)
  }
}

#' Calculate Rho Noise Projection and handle annihilation mode
#' @keywords internal
.ndx_calculate_rho_noise_projection <- function(workflow_state, current_pass_results, verbose) {
  
  opts_annihilation <- workflow_state$opts$opts_annihilation
  verbose_workflow <- workflow_state$opts$verbose
  
  # Build nuisance components list
  U_NDX_Nuisance_Combined_list <- list()
  if (!is.null(current_pass_results$rpca_components) && ncol(current_pass_results$rpca_components) > 0) {
    U_NDX_Nuisance_Combined_list$rpca <- current_pass_results$rpca_components
  }
  if (!is.null(current_pass_results$spectral_sines) && ncol(current_pass_results$spectral_sines) > 0) {
    U_NDX_Nuisance_Combined_list$spectral <- current_pass_results$spectral_sines
  }
  
  annihilation_active_this_pass <- FALSE
  
  # Handle annihilation mode orthogonalization
  if (opts_annihilation$annihilation_enable_mode && !is.null(workflow_state$U_GD_Lite_fixed_PCs) && 
      ncol(workflow_state$U_GD_Lite_fixed_PCs) > 0 && length(U_NDX_Nuisance_Combined_list) > 0) {
    
    U_NDX_Nuisance_Combined_list <- .ndx_orthogonalize_annihilation_components(
      U_NDX_Nuisance_Combined_list, workflow_state$U_GD_Lite_fixed_PCs, verbose_workflow
    )
    annihilation_active_this_pass <- TRUE
  }
  
  # Calculate rho noise projection
  rho_noise_projection <- .ndx_compute_rho_projection(U_NDX_Nuisance_Combined_list, workflow_state, verbose)
  
  list(
    rho_noise_projection = rho_noise_projection,
    annihilation_active_this_pass = annihilation_active_this_pass,
    U_NDX_Nuisance_Combined_list = U_NDX_Nuisance_Combined_list
  )
}

#' Orthogonalize components for annihilation mode
#' @keywords internal
.ndx_orthogonalize_annihilation_components <- function(U_NDX_Nuisance_Combined_list, U_GD_Lite_fixed_PCs, verbose_workflow) {
  
  if (verbose_workflow) message("  Annihilation Mode: Orthogonalizing NDX nuisance components against GLMdenoise-Lite PCs...")
  
  # Orthogonalize RPCA components
  if (!is.null(U_NDX_Nuisance_Combined_list$rpca) && ncol(U_NDX_Nuisance_Combined_list$rpca) > 0) {
    U_NDX_Nuisance_Combined_list$rpca_unique <- ndx_orthogonalize_matrix_against_basis(
      U_NDX_Nuisance_Combined_list$rpca, U_GD_Lite_fixed_PCs
    )
    U_NDX_Nuisance_Combined_list$rpca <- NULL
    if (verbose_workflow) {
      message(sprintf("    Orthogonalized %d RPCA components", ncol(U_NDX_Nuisance_Combined_list$rpca_unique)))
    }
  }
  
  # Orthogonalize spectral components
  if (!is.null(U_NDX_Nuisance_Combined_list$spectral) && ncol(U_NDX_Nuisance_Combined_list$spectral) > 0) {
    U_NDX_Nuisance_Combined_list$spectral_unique <- ndx_orthogonalize_matrix_against_basis(
      U_NDX_Nuisance_Combined_list$spectral, U_GD_Lite_fixed_PCs
    )
    U_NDX_Nuisance_Combined_list$spectral <- NULL
    if (verbose_workflow) {
      message(sprintf("    Orthogonalized %d spectral components", ncol(U_NDX_Nuisance_Combined_list$spectral_unique)))
    }
  }
  
  # Add GLMdenoise PCs
  U_NDX_Nuisance_Combined_list$gdlite <- U_GD_Lite_fixed_PCs
  if (verbose_workflow) {
    message(sprintf("    Added %d GLMdenoise-Lite PCs to nuisance set", ncol(U_GD_Lite_fixed_PCs)))
  }
  
  return(U_NDX_Nuisance_Combined_list)
}

#' Compute rho noise projection
#' @keywords internal
.ndx_compute_rho_projection <- function(U_NDX_Nuisance_Combined_list, workflow_state, verbose) {
  
  if (length(U_NDX_Nuisance_Combined_list) > 0 && !is.null(workflow_state$Y_residuals_current)) {
    U_NDX_Nuisance_Combined <- do.call(cbind, U_NDX_Nuisance_Combined_list)
    if (ncol(U_NDX_Nuisance_Combined) > 0) {
      P_noise_basis <- qr.Q(qr(U_NDX_Nuisance_Combined))
      Resid_proj_noise <- P_noise_basis %*% (crossprod(P_noise_basis, workflow_state$Y_residuals_current))
      var_resid_proj_noise <- sum(Resid_proj_noise^2, na.rm = TRUE)
      var_total_resid <- sum(workflow_state$Y_residuals_current^2, na.rm = TRUE)
      
      if (var_total_resid > 1e-9) {
        rho_value <- var_resid_proj_noise / var_total_resid
      } else {
        rho_value <- 0
      }
      
      if (verbose) message(sprintf("  Rho Noise Projection: %.4f", rho_value))
      return(rho_value)
    } else {
      if (verbose) message("  Rho Noise Projection: 0 (no combined nuisance components)")
      return(0)
    }
  } else {
    if (verbose) message("  Rho Noise Projection: NA (no nuisance components or residuals)")
    return(NA)
  }
}

#' Store additional diagnostics
#' @keywords internal
.ndx_store_additional_diagnostics <- function(pass_diagnostics, current_pass_results, annihilation_active_this_pass) {
  
  pass_diagnostics$k_rpca_global <- current_pass_results$k_rpca_global
  pass_diagnostics$num_hrf_clusters <- current_pass_results$num_hrf_clusters
  pass_diagnostics$num_spectral_sines <- current_pass_results$num_spectral_sines
  pass_diagnostics$lambda_parallel_noise <- current_pass_results$lambda_parallel_noise
  pass_diagnostics$lambda_perp_signal <- current_pass_results$lambda_perp_signal
  
  if (!is.null(current_pass_results$precision_weights)) {
    pass_diagnostics$precision_weight_summary <- stats::quantile(
      as.vector(current_pass_results$precision_weights),
      probs = c(0, 0.25, 0.5, 0.75, 1),
      na.rm = TRUE
    )
  } else {
    pass_diagnostics$precision_weight_summary <- NULL
  }
  
  pass_diagnostics$V_global_singular_values_from_rpca <- current_pass_results$V_global_singular_values_from_rpca
  
  pass_diagnostics$pass_options <- list(
    current_rpca_k_target = current_pass_results$k_rpca_global,
    annihilation_active_this_pass = annihilation_active_this_pass
  )
  
  return(pass_diagnostics)
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