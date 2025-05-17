#' Run a Single Pass of the ND-X Sprint 1 Denoising Workflow
#'
#' This function orchestrates the core modules of the ND-X pipeline developed in Sprint 1.
#' It performs initial residual generation, FIR HRF estimation, nuisance component
#' identification (RPCA and Spectral), AR(2) pre-whitening, and ridge regression.
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
#'   - `task_regressor_names_for_extraction` (character vector): Names of task regressors to extract betas for.
#' @param verbose Logical, if TRUE, print progress messages. Default: TRUE.
#'
#' @return A list containing key outputs from the workflow, such as:
#'   - `final_task_betas`: Extracted task betas after ridge regression.
#'   - `Y_residuals_whitened`: Residuals after whitening and ridge regression.
#'   - `Y_residuals_pass0`: Residuals from the initial GLM fit.
#'   - `ar_coeffs_voxelwise`: Voxelwise AR(2) coefficients.
#'   - `rpca_components`: Temporal components from RPCA.
#'   - `spectral_sines`: Sine/cosine components from spectral analysis.
#'   - `estimated_hrfs`: Table of estimated FIR HRFs.
#'   - `pass0_vars`: Baseline variance from initial GLM for DES calculation.
#'   - `na_mask_whitening`: Mask of NAs introduced by AR(2) whitening filter.
#'   - `X_full_design`: The constructed design matrix.
#' @importFrom fmrireg event_model design_matrix sampling_frame
#' @importFrom tibble is_tibble
#' @export
ndx_run_sprint1 <- function(Y_fmri,
                              events,
                              motion_params,
                              run_idx,
                              TR,
                              spike_TR_mask = NULL,
                              user_options = list(),
                              verbose = TRUE) {

  if (verbose) message("Starting ND-X Sprint 1 Workflow...")

  # --- 0. Validate inputs and merge default user_options --- 
  # (To be implemented: validate inputs, set up default options for each sub-module if not provided)
  
  # Default options for HRF estimation
  opts_pass0       <- user_options$opts_pass0 %||% list()

  default_opts_hrf <- list(
    hrf_fir_taps = 12L,
    hrf_fir_span_seconds = 24,
    good_voxel_R2_threshold = 0.05,
    cv_folds = 5L,
    lambda1_grid = 10^seq(-2, 1, length.out = 5),
    lambda2_grid = 10^seq(-3, 0, length.out = 5),
    hrf_min_good_voxels = 50L,
    return_full_model = FALSE
  )
  opts_hrf <- utils::modifyList(default_opts_hrf, user_options$opts_hrf %||% list())

  opts_rpca        <- user_options$opts_rpca %||% list()
  opts_spectral    <- user_options$opts_spectral %||% list()
  opts_whitening   <- user_options$opts_whitening %||% list()
  opts_ridge       <- user_options$opts_ridge %||% list()
  task_regressor_names <- user_options$task_regressor_names_for_extraction %||% character(0)
  
  # Initialize a list to store results
  results <- list()

  # --- 1. Initial GLM: Residual Generation (NDX-2) ---
  if (verbose) message("Step 1: Running Initial GLM for residual generation...")
  initial_glm_output <- ndx_initial_glm(
    Y_fmri = Y_fmri,
    events = events,
    motion_params = motion_params,
    run_idx = run_idx,
    TR = TR
  )
  results$Y_residuals_pass0 <- initial_glm_output$Y_residuals_current
  results$pass0_vars <- initial_glm_output$VAR_BASELINE_FOR_DES
  
  # --- 2. HRF Estimation (NDX-3) ---
  if (verbose) message("Step 2: Estimating initial FIR HRFs...")
  estimated_hrfs <- ndx_estimate_initial_hrfs(
    Y_fmri = Y_fmri, 
    pass0_residuals = initial_glm_output$Y_residuals_current, 
    events = events, 
    run_idx = run_idx, 
    TR = TR, 
    spike_TR_mask = spike_TR_mask, 
    user_options = opts_hrf
  )
  results$estimated_hrfs <- estimated_hrfs

  # --- 3. RPCA Nuisance Components (NDX-4) ---
  if (verbose) message("Step 3: Identifying RPCA nuisance components...")
  k_rpca_global <- opts_rpca$k_global_target %||% 5 # Example default
  rpca_components <- ndx_rpca_temporal_components_multirun(
    Y_residuals_cat = initial_glm_output$Y_residuals_current,
    run_idx = run_idx,
    k_global_target = k_rpca_global,
    user_options = opts_rpca
  )
  results$rpca_components <- rpca_components

  # --- 4. Spectral Nuisance Components (NDX-5) ---
  if (verbose) message("Step 4: Identifying spectral nuisance components...")
  if (!is.null(initial_glm_output$Y_residuals_current) && ncol(initial_glm_output$Y_residuals_current) > 0) {
    mean_residual_for_spectrum <- rowMeans(initial_glm_output$Y_residuals_current, na.rm = TRUE)
    spectral_sines <- ndx_spectral_sines(
      mean_residual_for_spectrum = mean_residual_for_spectrum,
      TR = TR,
      n_sine_candidates = opts_spectral$n_sine_candidates %||% 10,
      nyquist_guard_factor = opts_spectral$nyquist_guard_factor %||% 0.9,
      k_tapers = opts_spectral$k_tapers %||% 5,
      nw = opts_spectral$nw %||% 3
    )
    results$spectral_sines <- spectral_sines
  } else {
    if (verbose) message("  Skipping spectral analysis due to no Pass0 residuals or no voxels.")
    results$spectral_sines <- NULL
  }

  # --- 5. Construct Full Design Matrix for Whitening/Ridge --- 
  if (verbose) message("Step 5: Constructing full design matrix...")
  
  X_full_design <- ndx_build_design_matrix(
    estimated_hrfs = results$estimated_hrfs,
    events = events,
    motion_params = motion_params,
    rpca_components = results$rpca_components,
    spectral_sines = results$spectral_sines,
    run_idx = run_idx,
    TR = TR,
    poly_degree_val = opts_pass0$poly_degree, # poly_degree comes from opts_pass0
    verbose = verbose
  )
  results$X_full_design <- X_full_design
  
  # --- 6. AR(2) Pre-whitening (NDX-6) ---
  if (verbose) message("Step 6: Performing AR(2) pre-whitening...")
  if (!is.null(X_full_design)) {
    whitening_output <- ndx_ar2_whitening(
      Y_data = Y_fmri,
      X_design_full = X_full_design,
      Y_residuals_for_AR_fit = initial_glm_output$Y_residuals_current,
      order = opts_whitening$order %||% 2L,
      global_ar_on_design = opts_whitening$global_ar_on_design %||% TRUE
      # Removed deprecated max_ar_failures_prop
    )
    results$Y_whitened <- whitening_output$Y_whitened
    results$X_whitened <- whitening_output$X_whitened
    results$ar_coeffs_voxelwise <- whitening_output$ar_coeffs_voxelwise
    results$na_mask_whitening <- whitening_output$na_mask
  } else {
    if (verbose) message("  Skipping AR(2) whitening as X_full_design is NULL.")
    results$Y_whitened <- Y_fmri # Pass through unwhitened for ridge if X_full_design was missing
    results$X_whitened <- X_full_design # Will be NULL
    results$ar_coeffs_voxelwise <- NULL
    results$na_mask_whitening <- rep(FALSE, nrow(Y_fmri))
  }

  # --- 7. Ridge Regression (NDX-7) ---
  if (verbose) message("Step 7: Performing Ridge Regression...")
  if (!is.null(results$X_whitened) && !is.null(results$Y_whitened)) {
    ridge_betas_whitened <- ndx_solve_ridge(
      Y_whitened = results$Y_whitened,
      X_whitened = results$X_whitened,
      lambda_ridge = opts_ridge$lambda_ridge %||% 1.0, # Example default
      na_mask = results$na_mask_whitening
    )
    results$ridge_betas_whitened <- ridge_betas_whitened

    # Extract task betas
    if (!is.null(ridge_betas_whitened) && length(task_regressor_names) > 0) {
      final_task_betas <- ndx_extract_task_betas(
        betas_whitened = ridge_betas_whitened,
        X_whitened_colnames = colnames(results$X_whitened), 
        task_regressor_names = task_regressor_names,
        ar_coeffs_global = NULL # Unwhitening deferred beyond Sprint 1
      )
      results$final_task_betas <- final_task_betas
    } else {
      results$final_task_betas <- NULL
    }
    
    # Placeholder for whitened residuals after ridge
    # Y_residuals_whitened = Y_whitened - X_whitened %*% ridge_betas_whitened
    # results$Y_residuals_whitened <- Y_residuals_whitened 

  } else {
    if (verbose) message("  Skipping Ridge Regression due to missing whitened data/design.")
    results$ridge_betas_whitened <- NULL
    results$final_task_betas <- NULL
  }
  
  # --- 8. Finalize and Return --- 
  if (verbose) message("ND-X Sprint 1 Workflow finished.")
  return(results)
}

# Helper for default options, similar to base R's %||%
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
} 