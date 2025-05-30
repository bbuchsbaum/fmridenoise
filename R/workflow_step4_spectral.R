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
