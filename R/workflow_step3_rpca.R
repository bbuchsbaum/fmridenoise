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

