#' Process a Single Subject through the Full ND-X Iterative Denoising Workflow
#'
#' This function orchestrates the core modules of the ND-X pipeline. 
#' It iteratively performs initial residual generation, FIR HRF estimation, 
#' nuisance component identification (RPCA and Spectral), AR(2) pre-whitening, 
#' and ridge regression until convergence or max passes.
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

  # --- 0. Validate inputs and merge default user_options --- 
  # Default general workflow options
  max_passes <- user_options$max_passes %||% 3L
  min_des_gain_convergence <- user_options$min_des_gain_convergence %||% 0.005
  min_rho_noise_projection_convergence <- user_options$min_rho_noise_projection_convergence %||% 0.01

  # Default options for sub-modules (can be centralized or handled per module call)
  opts_pass0       <- user_options$opts_pass0 %||% list(poly_degree = 1) # Ensure poly_degree is available for build_design_matrix

  default_opts_hrf <- list(
    hrf_fir_taps = 12L,
    hrf_fir_span_seconds = 24,
    good_voxel_R2_threshold = 0.05,
    cv_folds = 5L,
    lambda1_grid = 10^seq(-2, 1, length.out = 5),
    lambda2_grid = 10^seq(-3, 0, length.out = 5),
    hrf_min_good_voxels = 50L,
    return_full_model = FALSE,
    hrf_cluster_method = "none", # Default to no clustering
    num_hrf_clusters = 1L,       # Consistent with no clustering
    hrf_cluster_merge_corr_thresh = 0.95,
    hrf_cluster_min_size = 5L,
    # Options for sparse event handling in HRF estimation (NDX-12)
    hrf_min_events_for_fir = 6L, 
    hrf_low_event_threshold = 12L,
    hrf_target_event_count_for_lambda_scaling = 20L,
    hrf_use_canonical_fallback_for_ultra_sparse = FALSE, # Default to attempting scaled/damped FIR
    # New options for project_hrf_cone
    hrf_cone_nonneg = TRUE,
    hrf_cone_unimodal = TRUE,
    hrf_cone_normalize_area = TRUE,
    verbose_hrf = FALSE # Added from previous refactor, ensure it's here
  )
  opts_hrf <- utils::modifyList(default_opts_hrf, user_options$opts_hrf %||% list())

  opts_rpca        <- user_options$opts_rpca %||% list()
  opts_spectral    <- user_options$opts_spectral %||% list()
  opts_whitening   <- user_options$opts_whitening %||% list()
  # opts_ridge       <- user_options$opts_ridge %||% list() # This was the old line
  # Default options for ridge, now more comprehensive for anisotropic and annihilation
  default_opts_ridge <- list(
    lambda_ridge = 1.0, # Isotropic lambda, for fallback or if anisotropic_ridge_enable=FALSE
    anisotropic_ridge_enable = TRUE,
    lambda_signal = 0.1, 
    lambda_noise_rpca = 1.0, 
    lambda_noise_spectral = 1.0,
    lambda_noise_other = 1.0,
    ridge_gcv_folds = 5L,
    ridge_update_lambda_aggressiveness = FALSE,
    # Lambdas for annihilation mode components - these are now part of the main ridge options
    lambda_noise_gdlite = 1.0, 
    lambda_noise_ndx_unique = 1.0 
  )
  opts_ridge <- utils::modifyList(default_opts_ridge, user_options$opts_ridge %||% list())

  # Default options for Annihilation mode (controls call to ndx_run_gdlite)
  default_opts_annihilation <- list(
    annihilation_enable_mode = FALSE,
    annihilation_gdlite_poly_degree = 1,
    annihilation_gdlite_k_max = 30L,
    annihilation_gdlite_r2_thresh_noise_pool = 0.05,
    annihilation_gdlite_tsnr_thresh_noise_pool = 30,
    annihilation_gdlite_r2_thresh_good_voxels = 0.05
  )
  opts_annihilation <- utils::modifyList(default_opts_annihilation, user_options$opts_annihilation %||% list())

  task_regressor_names <- user_options$task_regressor_names_for_extraction %||% character(0)
  verbose_workflow <- user_options$verbose %||% TRUE # Use a distinct name for workflow verbosity

  # Initialize storage for per-pass results
  diagnostics_per_pass <- list()
  beta_history_per_pass <- list()
  Y_residuals_current <- NULL # Will hold residuals from the previous pass
  VAR_BASELINE_FOR_DES <- NULL # From initial GLM
  current_overall_spike_TR_mask <- if (is.null(spike_TR_mask)) rep(FALSE, nrow(Y_fmri)) else as.logical(spike_TR_mask)
  # Placeholder for fixed GD-Lite PCs if Annihilation Mode is active
  U_GD_Lite_fixed_PCs <- NULL
  gdlite_initial_results <- NULL # To store full gdlite output if needed for diagnostics

  # --- Call ndx_run_gdlite if Annihilation Mode is enabled (once before the loop) ---
  if (opts_annihilation$annihilation_enable_mode) {
    if (verbose_workflow) message("Annihilation Mode enabled: Running GLMdenoise-Lite procedure upfront...")
    
    # Set verbose option for gdlite functions if main workflow verbose is on
    # This assumes gdlite functions use getOption("fmridenoise.verbose_gdlite", FALSE)
    original_gdlite_verbose_opt <- getOption("fmridenoise.verbose_gdlite")
    if (verbose_workflow) {
        options(fmridenoise.verbose_gdlite = TRUE)
    } else {
        options(fmridenoise.verbose_gdlite = FALSE)
    }

    gdlite_initial_results <- ndx_run_gdlite(
      Y_fmri = Y_fmri,
      events = events,
      run_idx = run_idx,
      TR = TR,
      motion_params = motion_params, # Pass original motion params
      poly_degree = opts_annihilation$annihilation_gdlite_poly_degree,
      k_max = opts_annihilation$annihilation_gdlite_k_max,
      r2_thresh_noise_pool = opts_annihilation$annihilation_gdlite_r2_thresh_noise_pool,
      tsnr_thresh_noise_pool = opts_annihilation$annihilation_gdlite_tsnr_thresh_noise_pool,
      r2_thresh_good_voxels = opts_annihilation$annihilation_gdlite_r2_thresh_good_voxels,
      perform_final_glm = FALSE # We only need the PCs for now
    )
    
    # Restore original gdlite verbose option
    options(fmridenoise.verbose_gdlite = original_gdlite_verbose_opt)

    if (!is.null(gdlite_initial_results) && !is.null(gdlite_initial_results$selected_pcs) && ncol(gdlite_initial_results$selected_pcs) > 0) {
      U_GD_Lite_fixed_PCs <- gdlite_initial_results$selected_pcs
      if (verbose_workflow) message(sprintf("  GLMdenoise-Lite selected %d PCs. These will be used as a fixed basis.", ncol(U_GD_Lite_fixed_PCs)))
    } else {
      if (verbose_workflow) message("  GLMdenoise-Lite did not select any PCs or failed. Annihilation mode will effectively be standard ND-X.")
      U_GD_Lite_fixed_PCs <- NULL # Ensure it's NULL if no PCs
    }
  } # End of upfront ndx_run_gdlite call

  # --- Iterative Refinement Loop (NDX-11) ---
  for (pass_num in 1:max_passes) {
    if (verbose) message(sprintf("\n--- Starting ND-X Pass %d/%d ---", pass_num, max_passes))
    
    current_pass_results <- list()

    # --- 1. Initial GLM (Pass 0 equivalent, or using residuals from previous ND-X pass) ---
    if (pass_num == 1) {
      if (verbose) message(sprintf("Pass %d: Running Initial GLM for residual generation...", pass_num))
      initial_glm_output <- ndx_initial_glm(
        Y_fmri = Y_fmri,
        events = events,
        motion_params = motion_params, # Initial GLM uses original motion params
        run_idx = run_idx,
        TR = TR
      )
      Y_residuals_current <- initial_glm_output$Y_residuals_current
      VAR_BASELINE_FOR_DES <- initial_glm_output$VAR_BASELINE_FOR_DES # Store for all DES calculations
      current_pass_results$Y_residuals_from_glm <- Y_residuals_current
      current_pass_results$pass0_vars <- VAR_BASELINE_FOR_DES
      pass0_residuals <- Y_residuals_current
    } else {
      if (verbose) message(sprintf("Pass %d: Using residuals from Pass %d.", pass_num, pass_num - 1))
      # Y_residuals_current is already set from the end of the previous pass
      if (is.null(Y_residuals_current)) {
          warning(sprintf("Pass %d: Y_residuals_current is NULL. Cannot proceed. Stopping iterations.", pass_num))
          break # Stop iteration
      }
    }
    
    # --- 2. HRF Estimation (NDX-3 / NDX-12) ---
    # For Sprint 2 NDX-12, this will become more complex with clustering and auto-adaptation.
    # For now, retain single global HRF estimation logic from Sprint 1.
    if (verbose) message(sprintf("Pass %d: Estimating FIR HRFs...", pass_num))
    estimated_hrfs_raw <- ndx_estimate_initial_hrfs(
      Y_fmri = Y_fmri,
      pass0_residuals = Y_residuals_current,
      events = events,
      run_idx = run_idx,
      TR = TR,
      spike_TR_mask = current_overall_spike_TR_mask,
      user_options = opts_hrf
    )
    current_pass_results$estimated_hrfs_per_cluster <- estimated_hrfs_raw # Store for diagnostics
    current_pass_results$num_hrf_clusters <-
      attr(estimated_hrfs_raw, "num_effective_clusters") %||% 1L

    if (!is.null(estimated_hrfs_raw) && tibble::is_tibble(estimated_hrfs_raw) && nrow(estimated_hrfs_raw) > 0) {
      if (verbose) message(sprintf("  Pass %d: Received %d raw HRF estimates (per cluster/condition).", pass_num, nrow(estimated_hrfs_raw)))
      # For ndx_build_design_matrix, select/derive one HRF per condition.
      # Strategy for Sprint 2 initial: Use HRF from cluster_id 1 if available, otherwise average or other rule.
      # Simplest for now: if num_hrf_clusters in opts_hrf is 1, it should be fine.
      # If > 1, take cluster_id 1, or average if that makes sense. 
      # For now, if num_hrf_clusters > 1, this will effectively use cluster 1 data for all voxels if no voxel-specific assignment is done before ndx_build_design_matrix.
      
      if (any(colnames(estimated_hrfs_raw) == "cluster_id")) {
         # If clustering was done and produced results for multiple clusters.
         # We need to select one HRF per *original condition name* for ndx_build_design_matrix.
         # Current strategy: pick cluster 1's HRF for each condition if available.
         estimated_hrfs_for_design <- estimated_hrfs_raw[estimated_hrfs_raw$cluster_id == 1L, 
                                                         c("condition", "hrf_estimate", "taps")]
         # Remove any potential duplicates if a condition somehow got multiple cluster_id=1 entries (unlikely)
         estimated_hrfs_for_design <- unique(tibble::as_tibble(estimated_hrfs_for_design))
      } else {
         # No cluster_id column, means it was global estimation (num_clusters=1)
         estimated_hrfs_for_design <- tibble::as_tibble(estimated_hrfs_raw[, c("condition", "hrf_estimate", "taps")])
      }
      
      original_conditions <- unique(as.character(events$condition))
      missing_hrfs <- original_conditions[!original_conditions %in% estimated_hrfs_for_design$condition]
      if (length(missing_hrfs) > 0 && verbose) {
          message(sprintf("    Pass %d: Not all conditions have a selected HRF for design matrix. Missing: %s", pass_num,
                          paste(missing_hrfs, collapse=", ")))
      }
      current_pass_results$estimated_hrfs <- estimated_hrfs_for_design
    } else {
      if (verbose) message(sprintf("Pass %d: HRF estimation returned NULL or empty. No task regressors will be built for design matrix.", pass_num))
      current_pass_results$estimated_hrfs <- NULL 
    }

    # --- 3. RPCA Nuisance Components (NDX-4 / NDX-13 / NDX-15) ---
    if (verbose) message(sprintf("Pass %d: Identifying RPCA nuisance components...", pass_num))

    # Prepare opts_rpca for the current pass, possibly with singular values from previous pass
    current_opts_rpca <- opts_rpca 
    if (pass_num > 1 && !is.null(diagnostics_per_pass[[pass_num - 1]]$V_global_singular_values_from_rpca)) {
      current_opts_rpca$adaptive_k_singular_values <- diagnostics_per_pass[[pass_num - 1]]$V_global_singular_values_from_rpca
      if (verbose) message(sprintf("    Pass %d: Using %d singular values from previous pass for RPCA rank adaptation.", 
                                 pass_num, length(current_opts_rpca$adaptive_k_singular_values)))
    } else if (pass_num > 1 && verbose) {
      message(sprintf("    Pass %d: No singular values from previous pass available for RPCA rank adaptation.", pass_num))
    }

    # Determine k_rpca_global (Adaptive rank logic, NDX-13)
    if (pass_num > 1 && !is.null(current_opts_rpca$adaptive_k_singular_values)) {
      sv_for_adapt <- current_opts_rpca$adaptive_k_singular_values
      if (is.numeric(sv_for_adapt) && length(sv_for_adapt) > 0) {
        drop_ratio <- current_opts_rpca$k_elbow_drop_ratio %||% 0.02
        k_min <- current_opts_rpca$k_rpca_min %||% 20L
        k_max <- current_opts_rpca$k_rpca_max %||% 50L
        # Ensure Auto_Adapt_RPCA_Rank is available (it's exported from ndx_rpca.R)
        k_rpca_global <- Auto_Adapt_RPCA_Rank(
          sv_for_adapt, 
          drop_ratio = drop_ratio,
          k_min = k_min,
          k_max = k_max
        )
        if (verbose) message(sprintf("    Pass %d: Adaptive k_rpca_global set to %d", pass_num, k_rpca_global))
      } else {
        k_rpca_global <- current_opts_rpca$k_global_target %||% 5 # Fallback
        if (verbose) message(sprintf("    Pass %d: Invalid/NULL singular values for adaptation, using k_rpca_global = %d from opts/default", pass_num, k_rpca_global))
      }
    } else {
      k_rpca_global <- current_opts_rpca$k_global_target %||% 5
      if (verbose && pass_num == 1) message(sprintf("    Pass %d: Using initial k_rpca_global = %d (from opts/default)", pass_num, k_rpca_global))
      else if (verbose && pass_num > 1) message(sprintf("    Pass %d: No adaptive singular values provided, using k_rpca_global = %d (from opts/default)", pass_num, k_rpca_global))
    }

    rpca_out <- ndx_rpca_temporal_components_multirun(
      Y_residuals_cat = Y_residuals_current,
      run_idx = run_idx,
      k_global_target = k_rpca_global,
      user_options = current_opts_rpca # Pass the potentially modified opts_rpca
    )
    current_pass_results$k_rpca_global <- k_rpca_global

    if (!is.null(rpca_out) && is.list(rpca_out)) {
      rpca_components <- rpca_out$C_components 
      if (!is.null(rpca_out$spike_TR_mask) && is.logical(rpca_out$spike_TR_mask) && 
          length(rpca_out$spike_TR_mask) == length(current_overall_spike_TR_mask)) {
        current_overall_spike_TR_mask <- current_overall_spike_TR_mask | rpca_out$spike_TR_mask 
      } else if (!is.null(rpca_out$spike_TR_mask) && verbose) {
        message("    Warning: spike_TR_mask from rpca_out is not a valid logical vector of correct length. Global spike mask not updated by this pass's RPCA.")
      }
      current_pass_results$S_matrix_rpca <- rpca_out$S_matrix_cat # Store S matrix
      current_pass_results$V_global_singular_values_from_rpca <- rpca_out$V_global_singular_values # Store for next pass

      if (!is.null(rpca_out$S_matrix_cat) &&
          is.matrix(rpca_out$S_matrix_cat) &&
          is.numeric(rpca_out$S_matrix_cat)) {
        current_pass_results$precision_weights <-
          ndx_precision_weights_from_S(rpca_out$S_matrix_cat)
        precision_weights_for_pass <- current_pass_results$precision_weights
      } else {
        if (verbose) message("    RPCA did not return a valid numeric S matrix; skipping precision weighting for this pass.")
        current_pass_results$precision_weights <- NULL
        precision_weights_for_pass <- NULL
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
      precision_weights_for_pass <- NULL
    }
    current_pass_results$rpca_components <- rpca_components
    current_pass_results$spike_TR_mask <- current_overall_spike_TR_mask 

    # --- 4. Spectral Nuisance Components (NDX-5 / NDX-14) ---
    # For Sprint 2 NDX-14, selection will be BIC/AIC based.
    if (verbose) message(sprintf("Pass %d: Identifying spectral nuisance components...", pass_num))
    if (!is.null(Y_residuals_current) && ncol(Y_residuals_current) > 0) {
      mean_residual_for_spectrum <- rowMeans(Y_residuals_current, na.rm = TRUE)
      spectral_sines <- ndx_spectral_sines(
        mean_residual_for_spectrum = mean_residual_for_spectrum,
        TR = TR,
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

    # --- 5. Construct Full Design Matrix (using ndx_build_design_matrix) --- 
    if (verbose) message(sprintf("Pass %d: Constructing full design matrix...", pass_num))
    
    # If in Annihilation Mode, need to handle the special orthogonalized components
    if (opts_annihilation$annihilation_enable_mode && !is.null(U_GD_Lite_fixed_PCs) && 
        ncol(U_GD_Lite_fixed_PCs) > 0) {
        
        # Standard components for base design
        X_full_design <- ndx_build_design_matrix(
          estimated_hrfs = current_pass_results$estimated_hrfs,
          events = events,
          motion_params = motion_params, # Original motion params are part of the base model always
          rpca_components = NULL,  # We'll add them separately with different prefixes
          spectral_sines = NULL,   # We'll add them separately with different prefixes
          run_idx = run_idx,
          TR = TR,
          poly_degree = opts_pass0$poly_degree, # From initial GLM setup
          verbose = FALSE # Reduce verbosity for internal calls
        )
        
        # Now handle the special components for Annihilation Mode
        if (!is.null(X_full_design)) {
            # Add GLMdenoise PCs with prefix 'gdlite_pc_'
            if (ncol(U_GD_Lite_fixed_PCs) > 0) {
                gdlite_pcols <- ncol(U_GD_Lite_fixed_PCs)
                colnames_gdlite <- paste0("gdlite_pc_", seq_len(gdlite_pcols))
                X_gdlite <- U_GD_Lite_fixed_PCs
                colnames(X_gdlite) <- colnames_gdlite
                X_full_design <- cbind(X_full_design, X_gdlite)
                if (verbose) message(sprintf("    Added %d GLMdenoise PCs to design matrix.", gdlite_pcols))
            }
            
            # Add orthogonalized RPCA components if they exist
            if (!is.null(U_NDX_Nuisance_Combined_list$rpca_unique) && 
                ncol(U_NDX_Nuisance_Combined_list$rpca_unique) > 0) {
                
                rpca_unique_pcols <- ncol(U_NDX_Nuisance_Combined_list$rpca_unique)
                colnames_rpca_unique <- paste0("rpca_unique_comp_", seq_len(rpca_unique_pcols))
                X_rpca_unique <- U_NDX_Nuisance_Combined_list$rpca_unique
                colnames(X_rpca_unique) <- colnames_rpca_unique
                X_full_design <- cbind(X_full_design, X_rpca_unique)
                if (verbose) message(sprintf("    Added %d orthogonalized RPCA components to design matrix.", rpca_unique_pcols))
            }
            
            # Add orthogonalized spectral components if they exist
            if (!is.null(U_NDX_Nuisance_Combined_list$spectral_unique) && 
                ncol(U_NDX_Nuisance_Combined_list$spectral_unique) > 0) {
                
                spectral_unique_pcols <- ncol(U_NDX_Nuisance_Combined_list$spectral_unique)
                colnames_spectral_unique <- paste0("spectral_unique_comp_", seq_len(spectral_unique_pcols))
                X_spectral_unique <- U_NDX_Nuisance_Combined_list$spectral_unique
                colnames(X_spectral_unique) <- colnames_spectral_unique
                X_full_design <- cbind(X_full_design, X_spectral_unique)
                if (verbose) message(sprintf("    Added %d orthogonalized spectral components to design matrix.", spectral_unique_pcols))
            }
        }
    } else {
        # Standard design matrix construction (non-Annihilation Mode)
        X_full_design <- ndx_build_design_matrix(
          estimated_hrfs = current_pass_results$estimated_hrfs,
          events = events,
          motion_params = motion_params, # Original motion params are part of the base model always
          rpca_components = current_pass_results$rpca_components, # From current pass
          spectral_sines = current_pass_results$spectral_sines,   # From current pass
          run_idx = run_idx,
          TR = TR,
          poly_degree = opts_pass0$poly_degree, # From initial GLM setup
          verbose = FALSE # Reduce verbosity for internal calls
        )
    }
    
    current_pass_results$X_full_design <- X_full_design
    
    # --- 6. AR(2) Pre-whitening (NDX-6) ---
    if (verbose) message(sprintf("Pass %d: Performing AR(2) pre-whitening...", pass_num))
    if (!is.null(X_full_design)) {
      # Residuals for AR fit should be from a GLM using X_full_design on Y_fmri (unwhitened)
      # This is a deviation from Sprint 1 where Y_residuals_pass0 was used.
      # For iterative refinement, AR model should be based on current best model residuals.
      temp_glm_for_ar_residuals <- tryCatch({
        stats::lm.fit(X_full_design, Y_fmri)$residuals
      }, error = function(e) {
        if(verbose) message(sprintf("  Pass %d: Error fitting temporary GLM for AR residuals: %s. Using Y_residuals_current for AR fit.", pass_num, e$message))
        Y_residuals_current # Fallback
      })
      
      whitening_output <- ndx_ar2_whitening(
        Y_data = Y_fmri,
        X_design_full = X_full_design,
        Y_residuals_for_AR_fit = temp_glm_for_ar_residuals,
        order = opts_whitening$order %||% 2L,
        global_ar_on_design = opts_whitening$global_ar_on_design %||% TRUE,
        weights = precision_weights_for_pass
      )
      current_pass_results$Y_whitened <- whitening_output$Y_whitened
      current_pass_results$X_whitened <- whitening_output$X_whitened
      current_pass_results$ar_coeffs_voxelwise <- whitening_output$ar_coeffs_voxelwise
      current_pass_results$na_mask_whitening <- whitening_output$na_mask
    } else {
      if (verbose) message(sprintf("  Pass %d: Skipping AR(2) whitening as X_full_design is NULL.", pass_num))
      current_pass_results$Y_whitened <- Y_fmri 
      current_pass_results$X_whitened <- X_full_design 
      current_pass_results$ar_coeffs_voxelwise <- NULL
      current_pass_results$na_mask_whitening <- rep(FALSE, nrow(Y_fmri))
    }

    # --- 7. Ridge Regression (NDX-7 / NDX-16 Anisotropic Ridge) ---
    if (verbose) message(sprintf("Pass %d: Performing Ridge/Anisotropic Regression...", pass_num))
    if (!is.null(current_pass_results$X_whitened) && !is.null(current_pass_results$Y_whitened) && ncol(current_pass_results$X_whitened) > 0) {
      
      n_regressors <- ncol(current_pass_results$X_whitened)

      if (opts_ridge$anisotropic_ridge_enable %||% TRUE) {
        regressor_names <- colnames(current_pass_results$X_whitened)
        if (!is.null(regressor_names)) {
          is_rpca_col <- grepl("^rpca_comp_", regressor_names)
          is_spectral_col <- grepl("^(sin_f|cos_f|spectral_comp_)", regressor_names)
          is_gdlite_col <- grepl("^gdlite_pc_", regressor_names)
          is_rpca_unique_col <- grepl("^rpca_unique_comp_", regressor_names)
          is_spectral_unique_col <- grepl("^spectral_unique_comp_", regressor_names)

          noise_rpca_col_indices <- which(is_rpca_col)
          noise_spectral_col_indices <- which(is_spectral_col)
          noise_gdlite_col_indices <- which(is_gdlite_col)
          noise_rpca_unique_col_indices <- which(is_rpca_unique_col)
          noise_spectral_unique_col_indices <- which(is_spectral_unique_col)

          I_reg <- diag(1, n_regressors)

          U_noise <- NULL
          if (length(noise_rpca_col_indices) > 0)
            U_noise <- cbind(U_noise, I_reg[, noise_rpca_col_indices, drop = FALSE])
          if (length(noise_spectral_col_indices) > 0)
            U_noise <- cbind(U_noise, I_reg[, noise_spectral_col_indices, drop = FALSE])

          U_gd <- NULL
          if (opts_annihilation$annihilation_enable_mode && length(noise_gdlite_col_indices) > 0)
            U_gd <- I_reg[, noise_gdlite_col_indices, drop = FALSE]

          U_unique <- NULL
          if (opts_annihilation$annihilation_enable_mode) {
            if (length(noise_rpca_unique_col_indices) > 0)
              U_unique <- cbind(U_unique, I_reg[, noise_rpca_unique_col_indices, drop = FALSE])
            if (length(noise_spectral_unique_col_indices) > 0)
              U_unique <- cbind(U_unique, I_reg[, noise_spectral_unique_col_indices, drop = FALSE])
          }

          proj_mats <- ndx_compute_projection_matrices(U_GD = U_gd, U_Unique = U_unique,
                                                       U_Noise = U_noise, n_regressors = n_regressors)

          res_var_est <- ndx_estimate_res_var_whitened(current_pass_results$Y_whitened,
                                                       current_pass_results$X_whitened,
                                                       current_pass_results$na_mask_whitening)

          lambda_ratio <- (opts_ridge$lambda_perp_signal %||% 0.1) /
                          (opts_ridge$lambda_parallel_noise %||% 10.0)
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

          current_pass_results$lambda_parallel_noise <- lambda_parallel_tuned
          current_pass_results$lambda_perp_signal <- lambda_ratio * lambda_parallel_tuned

          ridge_betas_whitened <- ndx_solve_anisotropic_ridge(
            Y_whitened = current_pass_results$Y_whitened,
            X_whitened = current_pass_results$X_whitened,
            projection_mats = proj_mats,
            lambda_values = lambda_values,
            na_mask = current_pass_results$na_mask_whitening,
            weights = precision_weights_for_pass,
            gcv_lambda = FALSE,
            res_var_scale = res_var_est
          )
        } else {
          ridge_betas_whitened <- NULL
          if (verbose) message("    Colnames for X_whitened not available, anisotropic ridge skipped.")
        }
      } else { # Isotropic ridge
        K_penalty_diag <- rep(opts_ridge$lambda_ridge %||% 1.0, n_regressors)
        current_pass_results$lambda_parallel_noise <- opts_ridge$lambda_ridge %||% 1.0
        current_pass_results$lambda_perp_signal <- opts_ridge$lambda_ridge %||% 1.0
        ridge_betas_whitened <- ndx_solve_anisotropic_ridge(
          Y_whitened = current_pass_results$Y_whitened,
          X_whitened = current_pass_results$X_whitened,
          K_penalty_diag = K_penalty_diag,
          na_mask = current_pass_results$na_mask_whitening,
          weights = precision_weights_for_pass
        )
      }
      current_pass_results$ridge_betas_whitened <- ridge_betas_whitened
      beta_history_per_pass[[pass_num]] <- ridge_betas_whitened 

      # Extract task betas
      if (!is.null(ridge_betas_whitened) && length(task_regressor_names) > 0) {
        final_task_betas_pass <- ndx_extract_task_betas(
          betas_whitened = ridge_betas_whitened,
          X_whitened_colnames = colnames(current_pass_results$X_whitened), 
          task_regressor_names = task_regressor_names,
          ar_coeffs_global = NULL 
        )
        current_pass_results$final_task_betas_pass <- final_task_betas_pass
      } else {
        current_pass_results$final_task_betas_pass <- NULL
      }
      
      # Calculate residuals for this pass (unwhitened, for DES and next pass input)
      if(!is.null(X_full_design) && !is.null(ridge_betas_whitened)) {
          Y_predicted_unwhitened = X_full_design %*% ridge_betas_whitened # Using original X_full_design
          Y_residuals_current = Y_fmri - Y_predicted_unwhitened # Update Y_residuals_current for next pass
          current_pass_results$Y_residuals_final_unwhitened <- Y_residuals_current
      } else {
          Y_residuals_current <- NULL # Cannot compute if parts are missing
          current_pass_results$Y_residuals_final_unwhitened <- Y_residuals_current
      }

    } else {
      if (verbose) message(sprintf("  Pass %d: Skipping Ridge Regression due to missing whitened data/design.", pass_num))
      current_pass_results$ridge_betas_whitened <- NULL
      current_pass_results$final_task_betas_pass <- NULL
      beta_history_per_pass[[pass_num]] <- NULL
      Y_residuals_current <- NULL # Reset if ridge failed
      current_pass_results$Y_residuals_final_unwhitened <- Y_residuals_current
    }
    
    # --- 8. Calculate Diagnostics for this Pass (NDX-9 / NDX-18) ---
    if (verbose) message(sprintf("Pass %d: Calculating diagnostics...", pass_num))
    pass_diagnostics <- list()

    # Ljung-Box whiteness test on residuals
    whitened_resids <- NULL
    if (!is.null(ridge_betas_whitened) &&
        !is.null(current_pass_results$Y_whitened) &&
        !is.null(current_pass_results$X_whitened)) {
      whitened_resids <- current_pass_results$Y_whitened -
        current_pass_results$X_whitened %*% ridge_betas_whitened
      if (!is.null(current_pass_results$na_mask_whitening)) {
        whitened_resids <- whitened_resids[!current_pass_results$na_mask_whitening, , drop = FALSE]
      }
    }
    pass_diagnostics$ljung_box_p <- ndx_ljung_box_pval(whitened_resids)
    
    # Calculate DES
    if (!is.null(Y_residuals_current) && !is.null(VAR_BASELINE_FOR_DES)) {
        pass_diagnostics$DES <- calculate_DES(current_residuals_unwhitened = Y_residuals_current, 
                                              VAR_BASELINE_FOR_DES = VAR_BASELINE_FOR_DES)
        if (verbose) message(sprintf("  Pass %d DES: %.4f", pass_num, pass_diagnostics$DES))
    } else {
        pass_diagnostics$DES <- NA
        if (verbose) message(sprintf("  Pass %d DES: NA (missing residuals or baseline variance)", pass_num))
    }
    
    # Calculate Rho Noise Projection
    U_NDX_Nuisance_Combined_list <- list()
    if (!is.null(current_pass_results$rpca_components) && ncol(current_pass_results$rpca_components) > 0) {
        U_NDX_Nuisance_Combined_list$rpca <- current_pass_results$rpca_components
    }
    if (!is.null(current_pass_results$spectral_sines) && ncol(current_pass_results$spectral_sines) > 0) {
        U_NDX_Nuisance_Combined_list$spectral <- current_pass_results$spectral_sines
    }
    
    # For Annihilation Mode: Orthogonalize NDX components against GLMdenoise PCs if available
    if (opts_annihilation$annihilation_enable_mode && !is.null(U_GD_Lite_fixed_PCs) && 
        ncol(U_GD_Lite_fixed_PCs) > 0 && length(U_NDX_Nuisance_Combined_list) > 0) {
        
        if (verbose_workflow) message("  Annihilation Mode: Orthogonalizing NDX nuisance components against GLMdenoise-Lite PCs...")
        
        # Orthogonalize each type of NDX component (RPCA, Spectral, etc.)
        if (!is.null(U_NDX_Nuisance_Combined_list$rpca) && ncol(U_NDX_Nuisance_Combined_list$rpca) > 0) {
            U_NDX_Nuisance_Combined_list$rpca_unique <- ndx_orthogonalize_matrix_against_basis(
                U_NDX_Nuisance_Combined_list$rpca, U_GD_Lite_fixed_PCs)
            U_NDX_Nuisance_Combined_list$rpca <- NULL # Remove original RPCA components
            if (verbose_workflow) message(sprintf("    Orthogonalized %d RPCA components", 
                                               ncol(U_NDX_Nuisance_Combined_list$rpca_unique)))
        }
        
        if (!is.null(U_NDX_Nuisance_Combined_list$spectral) && ncol(U_NDX_Nuisance_Combined_list$spectral) > 0) {
            U_NDX_Nuisance_Combined_list$spectral_unique <- ndx_orthogonalize_matrix_against_basis(
                U_NDX_Nuisance_Combined_list$spectral, U_GD_Lite_fixed_PCs)
            U_NDX_Nuisance_Combined_list$spectral <- NULL # Remove original spectral components
            if (verbose_workflow) message(sprintf("    Orthogonalized %d spectral components", 
                                               ncol(U_NDX_Nuisance_Combined_list$spectral_unique)))
        }
        
        # Add the fixed GLMdenoise PCs to the combined nuisance list
        U_NDX_Nuisance_Combined_list$gdlite <- U_GD_Lite_fixed_PCs
        if (verbose_workflow) message(sprintf("    Added %d GLMdenoise-Lite PCs to nuisance set", 
                                           ncol(U_GD_Lite_fixed_PCs)))
        
        # Update diagnostic to indicate annihilation was active this pass
        diagnostics_per_pass[[pass_num]]$pass_options$annihilation_active_this_pass <- TRUE
    }
    
    if (length(U_NDX_Nuisance_Combined_list) > 0 && !is.null(Y_residuals_current)) {
        U_NDX_Nuisance_Combined <- do.call(cbind, U_NDX_Nuisance_Combined_list)
        if (ncol(U_NDX_Nuisance_Combined) > 0) {
            P_noise_basis <- qr.Q(qr(U_NDX_Nuisance_Combined))
            Resid_proj_noise <- P_noise_basis %*% (crossprod(P_noise_basis, Y_residuals_current))
            var_resid_proj_noise <- sum(Resid_proj_noise^2, na.rm = TRUE)
            var_total_resid <- sum(Y_residuals_current^2, na.rm = TRUE)
            if (var_total_resid > 1e-9) { # Avoid division by zero if residuals are flat zero
                pass_diagnostics$rho_noise_projection <- var_resid_proj_noise / var_total_resid
            } else {
                pass_diagnostics$rho_noise_projection <- 0 # Or NA, if residuals are zero, no variance projects anywhere
            }
            if (verbose) message(sprintf("  Pass %d Rho Noise Projection: %.4f", pass_num, pass_diagnostics$rho_noise_projection %||% NA_real_))
        } else {
            pass_diagnostics$rho_noise_projection <- 0 # No nuisance components to project onto
            if (verbose) message(sprintf("  Pass %d Rho Noise Projection: 0 (no combined nuisance components)", pass_num))
        }
        
        # Store the latest nuisance components list for final output
        current_pass_results$U_NDX_Nuisance_Combined_list <- U_NDX_Nuisance_Combined_list
    } else {
        pass_diagnostics$rho_noise_projection <- NA # Cannot calculate if no nuisance or no residuals
        if (verbose) message(sprintf("  Pass %d Rho Noise Projection: NA (no nuisance components or residuals)", pass_num))
    }

    pass_diagnostics$k_rpca_global <- current_pass_results$k_rpca_global
    pass_diagnostics$num_hrf_clusters <- current_pass_results$num_hrf_clusters
    pass_diagnostics$num_spectral_sines <- current_pass_results$num_spectral_sines
    pass_diagnostics$lambda_parallel_noise <- current_pass_results$lambda_parallel_noise
    pass_diagnostics$lambda_perp_signal <- current_pass_results$lambda_perp_signal

    if (!is.null(precision_weights_for_pass)) {
        pass_diagnostics$precision_weight_summary <- stats::quantile(as.vector(precision_weights_for_pass),
                                                                     probs = c(0, 0.25, 0.5, 0.75, 1),
                                                                     na.rm = TRUE)
    } else {
        pass_diagnostics$precision_weight_summary <- NULL
    }
    
    pass_diagnostics$V_global_singular_values_from_rpca <- current_pass_results$V_global_singular_values_from_rpca # Add to pass diagnostics

    # Store options used for this pass if they can change adaptively or for reference
    diagnostics_per_pass[[pass_num]]$pass_options <- list(
      current_rpca_k_target = k_rpca_global,
      current_rpca_lambda = current_opts_rpca$rpca_lambda_fixed %||% current_opts_rpca$rpca_lambda_auto,
      # TODO: store spectral selection details, HRF clustering details etc.
      annihilation_active_this_pass = FALSE # Default, update if active
    )

    diagnostics_per_pass[[pass_num]] <- pass_diagnostics
    
    # --- 9. Convergence Check ---
    if (pass_num > 1) {
      converged_des <- FALSE
      converged_rho <- FALSE
      
      # DES gain check
      prev_DES <- diagnostics_per_pass[[pass_num - 1]]$DES
      current_DES <- pass_diagnostics$DES
      if (!is.na(prev_DES) && !is.na(current_DES)) {
        des_gain <- current_DES - prev_DES
        if (verbose) message(sprintf("  Pass %d DES gain: %.4f", pass_num, des_gain))
        if (des_gain < min_des_gain_convergence) {
          if (verbose) message(sprintf("  Convergence potentially met based on DES gain (%.4f < %.4f).", des_gain, min_des_gain_convergence))
          converged_des <- TRUE
        }
      } else {
        if (verbose) message("  Pass %d: Cannot check DES convergence due to NA DES values.")
        # If DES cannot be calculated, perhaps don't consider it converged based on DES.
      }
      
      # Rho noise projection check
      current_rho <- pass_diagnostics$rho_noise_projection
      if (!is.na(current_rho)) {
        if (verbose) message(sprintf("  Pass %d Current Rho Noise Projection: %.4f", pass_num, current_rho))
        if (current_rho < min_rho_noise_projection_convergence) {
          if (verbose) message(sprintf("  Convergence potentially met based on Rho Noise Projection (%.4f < %.4f).", current_rho, min_rho_noise_projection_convergence))
          converged_rho <- TRUE
        }
      } else {
          if (verbose) message("  Pass %d: Cannot check Rho convergence due to NA Rho value.")
      }
      
      # Stop if *either* convergence criterion is met and considered sufficient
      # Proposal: "Iteration stops if Denoising Efficacy Score (DES) gain is marginal (<0.5%) OR RHO?"
      # For now, let's assume if DES gain is too small, we stop. Rho is more for lambda adjustment.
      # Ticket NDX-11 says: "Implement convergence checks based on MIN_DES_GAIN_CONVERGENCE and MIN_RHO_NOISE_PROJECTION_CONVERGENCE"
      # This implies both could lead to stopping. Let's use OR for now: if DES gain is too small OR rho is small enough.
      if (converged_des || converged_rho) {
           if (verbose) message (sprintf("  Convergence criteria met (DES gain low: %s, Rho low: %s). Stopping iterations.", converged_des, converged_rho))
           break
      }

    }
    if (is.null(Y_residuals_current)) {
        if (verbose) message(sprintf("Pass %d: Y_residuals_current became NULL after ridge/residual calculation. Stopping iteration.", pass_num))
        break
    }

  } # End FOR loop (passes)
  
  # --- Finalize and Return --- 
  # Select results from the last successful pass or handle overall results structure
  final_results <- list(
     # Populate with outputs from the *last completed pass* stored in current_pass_results
     # or potentially the best pass based on DES if passes vary widely.
     # For now, just using last pass outputs.
     final_task_betas = current_pass_results$final_task_betas_pass,
     Y_residuals_final_unwhitened = current_pass_results$Y_residuals_final_unwhitened,
     ar_coeffs_voxelwise = current_pass_results$ar_coeffs_voxelwise,
     rpca_components = current_pass_results$rpca_components,
     spectral_sines = current_pass_results$spectral_sines,
     estimated_hrfs = current_pass_results$estimated_hrfs, # From last pass
     S_matrix_rpca_final = current_pass_results$S_matrix_rpca,
     pass0_vars = VAR_BASELINE_FOR_DES, # Original baseline variance
     pass0_residuals = pass0_residuals,
     na_mask_whitening = current_pass_results$na_mask_whitening,
     spike_TR_mask = current_overall_spike_TR_mask,
     X_full_design_final = current_pass_results$X_full_design, # From last pass
     diagnostics_per_pass = diagnostics_per_pass,
    beta_history_per_pass = beta_history_per_pass,
    num_passes_completed = pass_num
  )

  last_diag <- diagnostics_per_pass[[pass_num]]
  final_results$ljung_box_p <- last_diag$ljung_box_p
  final_results$num_hrf_clusters <- last_diag$num_hrf_clusters
  
  # Add Annihilation Mode specific results if active
  if (opts_annihilation$annihilation_enable_mode) {
    final_results$annihilation_mode_active <- TRUE
    final_results$gdlite_pcs <- U_GD_Lite_fixed_PCs
    
    # Add any orthogonalized components from the last pass
    if (!is.null(current_pass_results$U_NDX_Nuisance_Combined_list)) {
      final_results$rpca_orthogonalized <- current_pass_results$U_NDX_Nuisance_Combined_list$rpca_unique
      final_results$spectral_orthogonalized <- current_pass_results$U_NDX_Nuisance_Combined_list$spectral_unique
    }
    
    # Include GLMdenoise diagnostics
    if (!is.null(gdlite_initial_results)) {
      final_results$gdlite_diagnostics <- list(
        noise_pool_mask = gdlite_initial_results$noise_pool_mask,
        good_voxels_mask = gdlite_initial_results$good_voxels_mask,
        tsnr = gdlite_initial_results$tsnr,
        cv_r2 = gdlite_initial_results$cv_r2,
        optimal_k = gdlite_initial_results$optimal_k,
        r2_vals_by_k = gdlite_initial_results$r2_vals_by_k
      )
    }

    verdict_stats <- ndx_annihilation_verdict_stats(final_results)
    final_results$annihilation_var_ratio <- verdict_stats$var_ratio
    final_results$annihilation_verdict <- verdict_stats$verdict
  } else {
    final_results$annihilation_mode_active <- FALSE
  }
  
  if (verbose) message("ND-X Subject Processing finished.")
  return(final_results)
}

# Helper for default options, similar to base R's %||%
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
} 