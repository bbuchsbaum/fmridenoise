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

