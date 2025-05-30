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
