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
