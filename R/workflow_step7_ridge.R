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
