#' @keywords internal
#' @importFrom stats mad median quantile
# NO @export for project_hrf_cone (it's used by an exported function, but is internal itself)
project_hrf_cone <- function(h, 
                             nonneg = TRUE, 
                             unimodal = TRUE, 
                             normalize_area = TRUE,
                             verbose = FALSE) {
  
  if (!is.numeric(h) || length(h) == 0) {
    if (verbose) message("project_hrf_cone: input h is not numeric or is empty. Returning as is or empty.")
    return(h)
  }

  # NaN/NA Handling (Feedback 2.5)
  if (anyNA(h)) { 
    if(verbose) message("project_hrf_cone: Input h contains NA(s), setting to 0."); 
    h[is.na(h)] <- 0 
  }
  
  h_orig_for_debug <- h # Store after NA handling for all-zero check later
  original_unimodal_request <- unimodal # Store user's original unimodal request

  if (length(h) < 3L && unimodal) {
    if (verbose) message("project_hrf_cone: h length < 3, cannot enforce unimodality. Applying nonneg/norm if requested.")
    unimodal <- FALSE 
    if (original_unimodal_request && !unimodal) { # Warning if user asked for unimodal but it was skipped
        warning("project_hrf_cone: HRF length < 3, unimodal constraint skipped.", call. = FALSE)
    }
  }
  
  # 1. Force start at/near 0 by subtracting the first point
  if (length(h) > 0) { # Ensure h is not empty
    h <- h - h[1L]
  }
  
  # 2. Optional non-negativity constraint
  if (nonneg) {
    h <- pmax(h, 0)
  }
  
  # 3. Optional unimodality constraint
  pk_idx <- NULL # Define pk_idx outside to check its scope for the bump
  if (unimodal && length(h) >=3) {
    pk_idx <- which.max(h) 
    if (length(pk_idx) > 1) pk_idx <- pk_idx[1] 

    if (pk_idx > 1L) { # Ensure pk_idx is scalar and > 1 for first part
      h[1:pk_idx]  <- cummax(h[1:pk_idx])
    }
    if (pk_idx < length(h)) { # Ensure pk_idx is scalar and < length(h) for second part
      h[pk_idx:length(h)] <- rev(cummax(rev(h[pk_idx:length(h)])))
    }
    if (nonneg) h <- pmax(h, 0) # Re-apply nonneg after unimodal adjustments
  }

  # 4. Optional area normalization
  if (normalize_area) {
    s_h <- sum(h, na.rm = TRUE)
    if (abs(s_h) > 1e-9) { 
      h <- h / s_h
    } else {
      if (verbose && sum(abs(h_orig_for_debug), na.rm = TRUE) > 1e-9) {
        message("project_hrf_cone: HRF sum is near zero after projection, cannot normalize area to 1.")
      }
      # If sum is zero, and normalize_area is TRUE, what should h be? 
      # Current logic leaves h as is (likely all zeros or balanced positive/negative).
      # If nonneg=TRUE was also TRUE, then h should be all zeros here if s_h is zero.
    }
  }
  
  # Safeguard against all zeros if original was not all zeros
  current_sum_abs_h <- sum(abs(h), na.rm = TRUE)
  if (current_sum_abs_h < 1e-9 && sum(abs(h_orig_for_debug), na.rm = TRUE) > 1e-9 && length(h) > 0) {
      # If pk_idx wasn't set (e.g. unimodal=FALSE or length < 3), default peak_tap_idx
      peak_tap_idx <- if (!is.null(pk_idx) && length(pk_idx)==1) pk_idx else max(1L, floor(length(h) / 3L))
      if(peak_tap_idx > length(h) || peak_tap_idx < 1L) peak_tap_idx <- max(1L, floor(length(h)/3L))
      
      h[peak_tap_idx] <- 1e-6 
      if (normalize_area) { # If area was meant to be 1, and we just bumped it from zero
          # Re-normalize: h will be all zeros except for 1e-6 at peak_tap_idx, so sum is 1e-6
          # This results in h[peak_tap_idx] = 1, others 0.
          s_h_after_bump <- sum(h, na.rm = TRUE)
          if (abs(s_h_after_bump) > 1e-9) {
              h <- h / s_h_after_bump
          }
      }
      if(verbose) message(sprintf("project_hrf_cone: HRF became all zero, adding small bump at tap %d and re-evaluating normalization.", peak_tap_idx))
  }
  
  return(h)
}
