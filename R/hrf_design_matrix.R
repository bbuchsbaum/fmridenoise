#' @keywords internal
# NO @export for get_fir_design_matrix_for_condition
get_fir_design_matrix_for_condition <- function(condition_name, events_df,
                                                sampling_frame, fir_taps, TR,
                                                verbose = FALSE) {
  
  total_timepoints <- sum(sampling_frame$blocklens)
  
  # --- BEGIN DEBUG MESSAGES ---
  if (verbose) {
    message(sprintf("[get_fir_design_matrix_for_condition] Cond: %s", condition_name))
    message(sprintf("  Input fir_taps: %s (class: %s), Input TR: %s (class: %s)",
                    as.character(fir_taps), class(fir_taps), as.character(TR), class(TR)))
    message(sprintf("  sampling_frame$total_samples (total_timepoints): %s (class: %s)",
                    as.character(total_timepoints), class(total_timepoints)))
    if (!is.numeric(total_timepoints) || length(total_timepoints) != 1 || !is.finite(total_timepoints) || total_timepoints < 0) {
        message("  WARNING: total_timepoints is not a single positive finite number!")
    }
    if (!is.numeric(fir_taps) || length(fir_taps) != 1 || !is.finite(fir_taps) || fir_taps < 0) {
        message("  WARNING: fir_taps is not a single positive finite number!")
    }
  }
  # --- END DEBUG MESSAGES ---

  if (nrow(events_df) == 0 || fir_taps == 0) {
    # If no events or no FIR taps requested, return a zero matrix of correct dimensions
    # This is important for consistency if other conditions *do* produce regressors.
    return(matrix(0, nrow = total_timepoints, ncol = fir_taps))
  }
  
  # Manual FIR basis matrix construction (similar to .ndx_generate_task_regressors)
  X_fir_basis_cond <- matrix(0, nrow = total_timepoints, ncol = fir_taps)
  fir_colnames <- paste0(make.names(condition_name), "_FIRbasis", 0:(fir_taps-1)) # Unique basis names
  colnames(X_fir_basis_cond) <- fir_colnames

  current_run_offset <- 0 # To handle timepoints across concatenated runs
  for (run_idx_val in 1:length(sampling_frame$blocklens)) {
    run_length_tps <- sampling_frame$blocklens[run_idx_val]
    run_start_tp_global <- current_run_offset + 1
    run_end_tp_global <- current_run_offset + run_length_tps
    
    run_events <- events_df[events_df$blockids == run_idx_val, , drop = FALSE]

    if (nrow(run_events) > 0) {
      for (ev_idx in 1:nrow(run_events)) {
        event_onset_sec <- run_events$onsets[ev_idx]
        event_onset_tr_global_0idx <- floor(event_onset_sec / TR) + (run_start_tp_global - 1)

        for (tap_idx in 0:(fir_taps - 1)) { # 0-indexed tap
          target_timepoint_0idx <- event_onset_tr_global_0idx + tap_idx
          target_timepoint_1idx <- target_timepoint_0idx + 1

          if (target_timepoint_1idx >= run_start_tp_global && 
              target_timepoint_1idx <= run_end_tp_global && 
              target_timepoint_1idx <= total_timepoints) {
            X_fir_basis_cond[target_timepoint_1idx, tap_idx + 1] <- 1
          }
        }
      }
    }
    current_run_offset <- current_run_offset + run_length_tps
  }
  
  # The original function had checks for ncol(X_fir) != fir_taps
  # With manual construction, this should always match if fir_taps > 0.
  # If fir_taps was 0, we returned a 0-col matrix earlier.
  if (ncol(X_fir_basis_cond) != fir_taps) {
      # This case should ideally not be reached if fir_taps > 0 due to pre-allocation
      warning(sprintf("Manual FIR design for %s has %d columns, expected %d. This is unexpected.", 
                      condition_name, ncol(X_fir_basis_cond), fir_taps))
      # Fallback to a zero matrix of correct dimensions if something went very wrong
      return(matrix(0, nrow = total_timepoints, ncol = fir_taps, dimnames=list(NULL, fir_colnames))) 
  }
  
  return(X_fir_basis_cond)
}
