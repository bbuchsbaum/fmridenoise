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
