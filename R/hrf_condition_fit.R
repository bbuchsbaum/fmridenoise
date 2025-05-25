#' @keywords internal
# NO @export for estimate_hrf_for_condition
estimate_hrf_for_condition <- function(condition_name, events_for_condition,
                                       ybar_clean, block_ids_for_cv,
                                       overall_sampling_frame, valid_TRs_mask, TR,
                                       user_options) { 
  
  opts_hrf <- user_options 
  current_result <- list(condition = condition_name, hrf_estimate = NULL, taps = NULL, glmgen_fit = NULL)

  if (nrow(events_for_condition) == 0) {
      warning(paste("No events found for condition ", condition_name, " to estimate HRF.")) 
      return(current_result)
  }
  
  n_events_q <- nrow(events_for_condition) 
  lambda_scale_factor <- 1.0
  effective_hrf_fir_taps <- opts_hrf$hrf_fir_taps

  # Sparse event handling
  if (n_events_q < (opts_hrf$hrf_min_events_for_fir %||% 6L)) {
    if (opts_hrf$hrf_use_canonical_fallback_for_ultra_sparse %||% FALSE) {
      warning(sprintf("Condition '%s' has only %d events. Canonical HRF fallback not fully implemented; returning zero/damped FIR.", condition_name, n_events_q))
      current_result$hrf_estimate <- rep(0, opts_hrf$hrf_fir_taps)
      current_result$taps <- seq_len(opts_hrf$hrf_fir_taps)
      return(current_result)
    } else {
      warning(sprintf("Condition '%s' has only %d events (min_events_for_fir=%d). Returning zero/damped FIR.", 
                      condition_name, n_events_q, (opts_hrf$hrf_min_events_for_fir %||% 6L)))
      current_result$hrf_estimate <- rep(0, opts_hrf$hrf_fir_taps) 
      current_result$taps <- seq_len(opts_hrf$hrf_fir_taps)
      return(current_result)
    }
  } else if (n_events_q < (opts_hrf$hrf_low_event_threshold %||% 12L)) {
    target_events <- (opts_hrf$hrf_target_event_count_for_lambda_scaling %||% 20L)
    if (n_events_q > 0) { 
        lambda_scale_factor <- sqrt(target_events / n_events_q)
    }
    if (opts_hrf$verbose_hrf %||% FALSE) message(sprintf("  Condition '%s' has %d events. Applying lambda scaling factor: %.2f", 
                                 condition_name, n_events_q, lambda_scale_factor))
  }

  X_fir_cond_all_trs <- tryCatch({
       get_fir_design_matrix_for_condition(
          condition_name = condition_name,
          events_df = events_for_condition,
          sampling_frame = overall_sampling_frame,
          fir_taps = effective_hrf_fir_taps,
          TR = TR,
          verbose = (user_options$verbose_hrf %||% FALSE)
      )
  }, error = function(e) {
      warning(paste("Error generating FIR design for condition", condition_name, ":", e$message))
      return(NULL)
  })

  if (is.null(X_fir_cond_all_trs) || ncol(X_fir_cond_all_trs) == 0) {
    warning(paste("FIR design matrix for condition", condition_name, "is empty or NULL."))
    return(current_result)
  }
  
  X_fir_cond_clean <- X_fir_cond_all_trs[valid_TRs_mask, , drop = FALSE]
  
  if (nrow(X_fir_cond_clean) != length(ybar_clean)){
      stop(sprintf("Critical internal error: Mismatch in length of ybar_clean (%d) and rows of X_fir_cond_clean (%d) for condition: %s",
                 length(ybar_clean), nrow(X_fir_cond_clean), condition_name))
  }
  
  if (sum(abs(X_fir_cond_clean)) < 1e-6 ) { 
      warning(sprintf("No effective stimulus events for condition '%s' fall within non-spike TRs for FIR estimation.", condition_name))
      return(current_result)
  }

  min_trs_needed <- max(effective_hrf_fir_taps * 2, (opts_hrf$cv_folds %||% 2L) + 1) 
  if (nrow(X_fir_cond_clean) < min_trs_needed) { 
      warning(sprintf("Not enough valid TRs (%d, need ~%d) for FIR for condition: %s.", 
                      nrow(X_fir_cond_clean), min_trs_needed, condition_name))
      return(current_result)
  }
  
  cv_lambda1_grid <- (opts_hrf$lambda1_grid %||% 1) * lambda_scale_factor
  cv_lambda2_grid <- (opts_hrf$lambda2_grid %||% 0.1) * lambda_scale_factor 

  cv_results <- cv_fusedlasso(y = ybar_clean, X = X_fir_cond_clean,
                                lambda_grid = cv_lambda1_grid, 
                                gamma_grid = cv_lambda2_grid,  
                                k_folds = (opts_hrf$cv_folds %||% 2L),
                                block_ids = block_ids_for_cv)
  
  best_lambda1 <- cv_results$best_lambda 
  best_lambda2 <- cv_results$best_gamma  
  
  if (opts_hrf$verbose_hrf %||% FALSE) {
    message(sprintf("  Condition '%s': Best lambda1 (L1 diff): %.4f, Best lambda2 (L2 coeff): %.4f (from scaled grids if applicable, scale factor: %.2f)", 
                  condition_name, best_lambda1, best_lambda2, lambda_scale_factor))
  }
  
  X_mat <- X_fir_cond_clean
  p_coeffs_final <- ncol(X_mat)
  D_1d_final <- NULL
  if (p_coeffs_final > 1) {
    D_1d_final <- diff(diag(p_coeffs_final), differences = 1)
  } else if (p_coeffs_final == 1) {
    D_1d_final <- matrix(0, nrow=0, ncol=1) 
  }
  
  final_fit_path <- tryCatch({
    genlasso::fusedlasso(y = ybar_clean, X = X_mat, D = D_1d_final, gamma = best_lambda2)
  }, error = function(e) {
    warning(sprintf("Final genlasso::fusedlasso fit failed for condition '%s' with gamma=%.4f: %s", 
                    condition_name, best_lambda2, e$message))
    return(NULL)
  })
  
  if (is.null(final_fit_path) || !inherits(final_fit_path, "genlasso")) {
    warning(sprintf("Final genlasso::fusedlasso fit for condition '%s' is NULL or not a genlasso object.", condition_name))
    return(current_result)
  }
  
  fir_coeffs <- tryCatch({
    stats::coef(final_fit_path, lambda = best_lambda1)$beta
  }, error = function(e) {
    warning(sprintf("Failed to extract coefficients for condition '%s' with lambda=%.4f from genlasso fit: %s", 
                    condition_name, best_lambda1, e$message))
    return(NULL)
  })
  
  if (is.null(fir_coeffs)) {
    return(current_result) 
  }
  
  if (length(fir_coeffs) != effective_hrf_fir_taps) {
    if (opts_hrf$verbose_hrf %||% FALSE) {
        message(sprintf("    Condition %s: Coeff length %d from fusedlasso, expected %d. Padding/truncating.",
                        condition_name, length(fir_coeffs), effective_hrf_fir_taps)) 
    }
    temp_coeffs <- rep(0, effective_hrf_fir_taps)
    len_to_copy <- min(length(fir_coeffs), effective_hrf_fir_taps)
    if (len_to_copy > 0) temp_coeffs[1:len_to_copy] <- fir_coeffs[1:len_to_copy]
    fir_coeffs <- temp_coeffs
  }
  
  # Apply the new cone projection
  current_result$hrf_estimate <- project_hrf_cone(
                                      fir_coeffs, 
                                      nonneg = opts_hrf$hrf_cone_nonneg %||% TRUE, 
                                      unimodal = opts_hrf$hrf_cone_unimodal %||% TRUE, 
                                      normalize_area = opts_hrf$hrf_cone_normalize_area %||% TRUE,
                                      verbose = (opts_hrf$verbose_hrf %||% FALSE)
                                    )
  current_result$taps <- seq_len(effective_hrf_fir_taps)
  if (opts_hrf$return_full_model) {
      current_result$glmgen_fit <- final_fit_path
  }
  return(current_result)
}
