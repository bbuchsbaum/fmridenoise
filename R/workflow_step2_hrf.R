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
