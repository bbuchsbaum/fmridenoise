#' Identify Sinusoidal Nuisance Regressors from a Spectrum
#'
#' This function uses multi-taper spectral estimation to identify prominent
#' sinusoidal components in a residual time series (typically a mean residual).
#' It then generates pairs of sine and cosine regressors at the frequencies
#' of the identified spectral peaks.
#'
#' @param mean_residual_for_spectrum A numeric vector representing the time series
#'   (e.g., mean residuals across voxels) from which to estimate the spectrum.
#' @param TR Numeric, the repetition time of the fMRI data in seconds.
#' @param n_sine_candidates Integer, the maximum number of peak frequencies to
#'   convert into sine/cosine regressors. Defaults to 6.
#' @param nyquist_guard_factor Numeric, a factor (0 to 1) to limit peak searching
#'   to frequencies below this fraction of the Nyquist frequency. Helps avoid
#'   aliasing or unstable estimates near Nyquist. Defaults to 0.9.
#' @param k_tapers Integer, the number of tapers to use in `multitaper::spec.mtm`.
#'   Defaults to 5.
#' @param nw Numeric, the time-bandwidth product for `multitaper::spec.mtm`.
#'   Defaults to 3.
#' @param selection_criterion Character string, either "BIC" or "AIC". Default "BIC".
#' @param selection_delta_threshold Numeric, minimum decrease in the chosen information
#'   criterion required to retain a candidate pair. Default 2.
#' @param verbose Logical, whether to print verbose diagnostic messages.
#'
#' @return A matrix with `2 * n_selected_peaks` columns and
#'   `length(mean_residual_for_spectrum)` rows. Each pair of columns
#'   represents sine and cosine regressors for an identified peak
#'   frequency. Attributes \code{"freq_hz"} and \code{"freq_rad_s"}
#'   store the selected peak frequencies in Hz and radians/sec,
#'   respectively. If no peaks survive selection, a zero-column matrix
#'   with these attributes is returned. The function returns \code{NULL}
#'   only when the input is invalid or when spectral estimation fails.
#'
#' @examples
#' \dontrun{
#'   set.seed(42)
#'   TR_val <- 2.0
#'   n_timepoints <- 200
#'   time_points <- seq(0, (n_timepoints - 1) * TR_val, by = TR_val)
#'
#'   # Create a synthetic residual with a 0.05 Hz sine wave + noise
#'   freq_of_interest <- 0.05 # Hz
#'   synthetic_signal <- sin(2 * pi * freq_of_interest * time_points)
#'   noise <- rnorm(n_timepoints, sd = 0.5)
#'   mean_resid <- synthetic_signal + noise
#'
#'   U_sines <- ndx_spectral_sines(mean_resid, TR = TR_val, n_sine_candidates = 3)
#'   if (!is.null(U_sines)) {
#'     print(paste("Generated", ncol(U_sines) / 2, "sine/cosine pairs."))
#'     print("Identified frequencies (Hz):")
#'     print(attr(U_sines, "freq_hz"))
#'     # plot(time_points, U_sines[,1], type='l', main="First Sine Regressor")
#'   }
#' }
#'
#' @importFrom multitaper spec.mtm
#' @importFrom pracma findpeaks
#' @export
ndx_spectral_sines <- function(mean_residual_for_spectrum, TR,
                               n_sine_candidates = 6,
                               nyquist_guard_factor = 0.9,
                               k_tapers = 5, nw = 3,
                               selection_criterion = "BIC",
                               selection_delta_threshold = 2,
                               verbose = FALSE) {
  val <- .validate_spectral_inputs(mean_residual_for_spectrum, TR,
                                   n_sine_candidates,
                                   selection_delta_threshold)
  if (!is.null(val$error)) return(val$error)
  mean_residual_for_spectrum <- val$vec
  n_sine_candidates <- val$n_sine_candidates
  selection_delta_threshold <- val$selection_delta_threshold

  adj <- .adjust_mtm_params(k_tapers, nw,
                            length(mean_residual_for_spectrum), verbose)
  if (!is.null(adj$error)) return(adj$error)
  k_tapers <- adj$k_tapers
  nw <- adj$nw

  mt_res <- .estimate_spectrum(mean_residual_for_spectrum, TR, k_tapers, nw)
  if (is.null(mt_res)) {
    return(.empty_spec_matrix(length(mean_residual_for_spectrum)))
  }

  freqs_hz <- .find_spectral_peaks(mt_res, TR, nyquist_guard_factor,
                                   n_sine_candidates, verbose)
  if (length(freqs_hz) == 0) {
    if (verbose) message("No peaks selected after uniqueness filter to generate sine/cosine pairs.")
    return(.empty_spec_matrix(length(mean_residual_for_spectrum)))
  }

  U_spectral_sines <- .build_sine_regressors(freqs_hz,
                                             length(mean_residual_for_spectrum),
                                             TR)

  if (verbose) {
    message(sprintf("Generated %d sine/cosine pairs from %d spectral peaks.",
                    ncol(U_spectral_sines) / 2, length(freqs_hz)))
  }

  selected_regressors <- Select_Significant_Spectral_Regressors(
    y = mean_residual_for_spectrum,
    U_candidates = U_spectral_sines,
    criterion = selection_criterion,
    delta_threshold = selection_delta_threshold,
    verbose = verbose
  )

  return(selected_regressors)
}

#' Create an empty matrix for spectral regressors
#'
#' @keywords internal
.empty_spec_matrix <- function(n_rows) {
  res <- matrix(numeric(0), nrow = n_rows, ncol = 0,
                dimnames = list(NULL, character(0)))
  attr(res, "freq_hz") <- numeric(0)
  attr(res, "freq_rad_s") <- numeric(0)
  attr(res, "selected_pair_indices") <- integer(0)
  res
}

#' Validate and preprocess inputs for ndx_spectral_sines
#'
#' Returns a list containing a centered time series, the (possibly adjusted)
#' `n_sine_candidates`, `selection_delta_threshold`, and an `error` element. If
#' validation fails the `error` element contains the object that should be
#' returned by `ndx_spectral_sines`.
#'
#' @keywords internal
.validate_spectral_inputs <- function(mean_residual_for_spectrum, TR,
                                      n_sine_candidates,
                                      selection_delta_threshold) {
  if (!is.numeric(mean_residual_for_spectrum) ||
      length(mean_residual_for_spectrum) == 0) {
    warning("mean_residual_for_spectrum must be a non-empty numeric vector.")
    return(list(error = NULL))
  }
  if (!is.numeric(TR) || TR <= 0) {
    warning("TR must be a positive numeric value.")
    return(list(error = NULL))
  }

  if (!all(is.finite(mean_residual_for_spectrum))) {
    warning("mean_residual_for_spectrum contains non-finite values.")
    return(list(error = .empty_spec_matrix(length(mean_residual_for_spectrum))))
  }

  default_n_sine_candidates <- 6
  if (!is.numeric(n_sine_candidates) || length(n_sine_candidates) != 1 ||
      n_sine_candidates <= 0) {
    warning(sprintf("n_sine_candidates must be > 0. Using default %d.",
                    default_n_sine_candidates))
    n_sine_candidates <- default_n_sine_candidates
  }

  default_selection_delta_threshold <- 2
  if (!is.numeric(selection_delta_threshold) ||
      length(selection_delta_threshold) != 1 ||
      selection_delta_threshold < 0) {
    warning(sprintf(
      "selection_delta_threshold must be >= 0. Using default %.2f.",
      default_selection_delta_threshold
    ))
    selection_delta_threshold <- default_selection_delta_threshold
  }

  vec <- base::scale(mean_residual_for_spectrum, center = TRUE,
                     scale = FALSE)[, 1]

  list(vec = vec,
       n_sine_candidates = n_sine_candidates,
       selection_delta_threshold = selection_delta_threshold,
       error = NULL)
}

#' Adjust multitaper parameters for series length
#'
#' Ensures `nw` and `k_tapers` are compatible with the data length. If the
#' resulting values are invalid a list with an `error` element is returned.
#'
#' @keywords internal
.adjust_mtm_params <- function(k_tapers, nw, series_length, verbose = FALSE) {
  original_nw <- nw
  nw <- min(nw, floor((series_length - 1) / 2))
  if (nw < original_nw && verbose) {
    message(sprintf(
      "  ndx_spectral_sines: nw was reduced from %s to %s to be compatible with series length %d.",
      original_nw, nw, series_length))
  }
  if (nw <= 0) {
    warning(sprintf(
      "nw became %.2f after adjustment. Must be > 0. Aborting spectral step.",
      nw))
    return(list(error = .empty_spec_matrix(series_length)))
  }

  original_k_tapers <- k_tapers
  k_tapers <- as.integer(min(k_tapers, 2 * nw - 1, series_length - 1))
  if (k_tapers < original_k_tapers && verbose) {
    message(sprintf(
      "  ndx_spectral_sines: k_tapers reduced from %d to %d based on nw (%.2f) and series length (%d).",
      original_k_tapers, k_tapers, nw, series_length))
  }

  if (k_tapers < 1) {
    warning(sprintf(
      "k_tapers became %d after adjustment (nw=%.2f, series_length=%d). Must be >= 1. Aborting spectral step.",
      k_tapers, nw, series_length))
    return(list(error = .empty_spec_matrix(series_length)))
  }

  list(k_tapers = k_tapers, nw = nw, error = NULL)
}

#' Run multitaper spectrum estimation
#'
#' Wrapper around `multitaper::spec.mtm`. Returns `NULL` when estimation fails.
#'
#' @keywords internal
.estimate_spectrum <- function(ts, TR, k_tapers, nw) {
  mt_res <- NULL
  tryCatch({
    mt_res <- multitaper::spec.mtm(ts,
                                   k = k_tapers,
                                   nw = nw,
                                   deltat = TR,
                                   jackknife = FALSE,
                                   returnInternals = FALSE,
                                   plot = FALSE)
  }, error = function(e) {
    warning(paste("multitaper::spec.mtm failed:", e$message))
    assign("mt_res", NULL, envir = parent.env(environment()))
  })

  if (is.null(mt_res) || is.null(mt_res$spec) || is.null(mt_res$freq) ||
      length(mt_res$spec) == 0) {
    warning("Spectrum estimation via spec.mtm did not yield valid spec or freq.")
    return(NULL)
  }
  mt_res
}

#' Identify prominent spectral peaks
#'
#' Searches for peaks below a guarded Nyquist frequency and returns the selected
#' frequencies in Hz.
#'
#' @keywords internal
.find_spectral_peaks <- function(mt_res, TR, nyquist_guard_factor,
                                 n_sine_candidates, verbose = FALSE) {
  nyquist_freq <- 1 / (2 * TR)
  guard_factor <- max(0, min(1, nyquist_guard_factor))

  freq_upper_bound <- guard_factor * nyquist_freq
  keep_indices <- mt_res$freq > 0 & mt_res$freq <= freq_upper_bound

  if (sum(keep_indices) == 0) {
    warning("No frequencies to search for peaks after applying guard factor and excluding DC.")
    return(numeric(0))
  }

  spec_to_search <- mt_res$spec[keep_indices]
  freq_to_search <- mt_res$freq[keep_indices]

  if (verbose) {
    message(sprintf(
      "  Spectral search range: %.4f Hz to %.4f Hz (%d points)",
      min(freq_to_search), max(freq_to_search), length(freq_to_search)))
    message(sprintf(
      "  Summary of spec_to_search: Min=%.2e, Med=%.2e, Mean=%.2e, Max=%.2e, SD=%.2e",
      min(spec_to_search), median(spec_to_search), mean(spec_to_search),
      max(spec_to_search), sd(spec_to_search)))
  }

  peak_info_all <- pracma::findpeaks(spec_to_search, nups = 1, ndowns = 1,
                                     sortstr = TRUE)
  if (verbose && !is.null(peak_info_all)) {
    message(sprintf(
      "  findpeaks found %d initial peaks. Top peak amplitude: %.2e at index %d",
      nrow(peak_info_all),
      if (nrow(peak_info_all) > 0) peak_info_all[1, 1] else NA,
      if (nrow(peak_info_all) > 0) peak_info_all[1, 2] else NA))
  } else if (verbose && is.null(peak_info_all)) {
    message("  findpeaks returned NULL.")
  }

  if (is.null(peak_info_all) || nrow(peak_info_all) == 0) {
    if (verbose) message("No initial peaks found in the spectrum by pracma::findpeaks.")
    return(numeric(0))
  }

  prom_thresh <- median(spec_to_search, na.rm = TRUE) +
    3 * mad(spec_to_search, na.rm = TRUE)
  if (verbose) message(sprintf("  Prominence threshold for peaks: %.2e", prom_thresh))

  if (!is.finite(prom_thresh)) {
    prom_thresh <- -Inf
    warning("Prominence threshold for peak picking was not finite (e.g. MAD was zero). Filter effectively disabled.")
  }

  peak_info <- peak_info_all[peak_info_all[, 1] > prom_thresh, , drop = FALSE]
  if (verbose && !is.null(peak_info)) {
    message(sprintf(
      "  After prominence filter, %d peaks remain. Top peak amplitude: %.2e",
      nrow(peak_info), if (nrow(peak_info) > 0) peak_info[1, 1] else NA))
  } else if (verbose && is.null(peak_info)) {
    message("  peak_info is NULL after prominence filter (should be 0-row matrix if no peaks pass).")
  }

  if (is.null(peak_info) || nrow(peak_info) == 0) {
    if (verbose) message(sprintf(
      "No significant peaks found after prominence filter (threshold: %.4g).",
      prom_thresh))
    return(numeric(0))
  }

  num_peaks_to_select <- min(n_sine_candidates, nrow(peak_info))
  sel_idx <- peak_info[1:num_peaks_to_select, 2]

  selected_frequencies_hz <- freq_to_search[sel_idx]
  if (length(selected_frequencies_hz) > 0) {
    orig_n <- length(selected_frequencies_hz)
    selected_frequencies_hz <- unique(selected_frequencies_hz)
    if (length(selected_frequencies_hz) < orig_n && verbose) {
      message(sprintf(
        "  Reduced to %d unique peak frequencies from %d for regressor generation.",
        length(selected_frequencies_hz), orig_n))
    }
  }

  selected_frequencies_hz
}

#' Build sine and cosine regressors for the given frequencies
#'
#' Returns a matrix with sine columns followed by cosine columns. Attributes
#' `freq_hz` and `freq_rad_s` store the input frequencies.
#'
#' @keywords internal
.build_sine_regressors <- function(freqs_hz, series_length, TR) {
  if (length(freqs_hz) == 0) {
    return(.empty_spec_matrix(series_length))
  }

  omega_rad_s <- 2 * pi * freqs_hz
  time_vector_sec <- (seq_len(series_length) - 1) * TR

  U_sin_raw <- sapply(omega_rad_s, function(w) sin(w * time_vector_sec))
  U_cos_raw <- sapply(omega_rad_s, function(w) cos(w * time_vector_sec))

  if (is.vector(U_sin_raw)) U_sin_raw <- matrix(U_sin_raw, ncol = 1)
  if (is.vector(U_cos_raw)) U_cos_raw <- matrix(U_cos_raw, ncol = 1)

  U_sin <- apply(U_sin_raw, 2, function(col) {
    nrm <- sqrt(sum(col^2))
    if (nrm > .Machine$double.eps) col / nrm else col
  })
  U_cos <- apply(U_cos_raw, 2, function(col) {
    nrm <- sqrt(sum(col^2))
    if (nrm > .Machine$double.eps) col / nrm else col
  })

  if (is.vector(U_sin)) U_sin <- matrix(U_sin, ncol = 1)
  if (is.vector(U_cos)) U_cos <- matrix(U_cos, ncol = 1)

  colnames(U_sin) <- paste0("sin_f", sprintf("%.4f", freqs_hz))
  colnames(U_cos) <- paste0("cos_f", sprintf("%.4f", freqs_hz))

  U <- cbind(U_sin, U_cos)
  attr(U, "freq_hz") <- freqs_hz
  attr(U, "freq_rad_s") <- omega_rad_s
  U
}

#' Select Spectral Regressors via Information Criterion
#'
#' Given a set of candidate sine/cosine regressors, iteratively add the pair
#' that yields the largest improvement in BIC (or AIC) when regressing the
#' provided time-series. Addition stops when no remaining pair improves the
#' criterion by more than `delta_threshold`.
#'
#' @param y Numeric vector of data to be modeled (typically the mean residual
#'   time series used to estimate the spectrum).
#' @param U_candidates Matrix of sine/cosine regressors. Columns should come in
#'   pairs (sine then cosine) corresponding to frequencies.
#' @param criterion Character string, either "BIC" or "AIC". Default "BIC".
#' @param delta_threshold Numeric, minimum decrease in the chosen information
#'   criterion required to retain a candidate pair. Default 2.
#' @param verbose Logical, whether to print verbose diagnostic messages.
#' @return Matrix of selected regressors (or NULL if none selected). Attributes
#'   "freq_hz" and "freq_rad_s" are propagated for the selected frequencies.
#' @keywords internal
#' @export
Select_Significant_Spectral_Regressors <- function(y, U_candidates,
                                                   criterion = c("BIC", "AIC"),
                                                   delta_threshold = 2,
                                                   verbose = FALSE) {
  if (is.null(U_candidates) || !is.matrix(U_candidates)) { # Allow 0-col matrix from caller
    if (verbose) message("[Select_SSR] U_candidates is NULL or not a matrix. Returning empty attributed matrix.")
    return(.empty_spec_matrix(length(y)))
  }
  if (ncol(U_candidates) < 2 && ncol(U_candidates) != 0 ) { # if not 0, must be pairs
     warning("[Select_SSR] U_candidates columns not in pairs or empty. Returning empty matrix.")
    return(.empty_spec_matrix(length(y)))
  }
  
  criterion <- match.arg(criterion)

  default_delta_threshold <- 2
  if (!is.numeric(delta_threshold) || length(delta_threshold) != 1 ||
      delta_threshold < 0) {
    warning(sprintf("[Select_SSR] delta_threshold must be >= 0. Using default %.2f.",
                    default_delta_threshold))
    delta_threshold <- default_delta_threshold
  }

  n_pairs <- ncol(U_candidates) %/% 2
  
  # If U_candidates is a 0-column matrix (e.g. no peaks from ndx_spectral_sines), n_pairs will be 0
  if (n_pairs < 1) {
    if (verbose) message("[Select_SSR] No candidate pairs to select from. Returning empty attributed matrix.")
    res_empty <- .empty_spec_matrix(length(y))
    hz_attr <- attr(U_candidates, "freq_hz")
    rad_attr <- attr(U_candidates, "freq_rad_s")
    attr(res_empty, "freq_hz") <- if (!is.null(hz_attr)) hz_attr else numeric(0)
    attr(res_empty, "freq_rad_s") <- if (!is.null(rad_attr)) rad_attr else numeric(0)
    return(res_empty)
  }

  calc_ic <- function(y, X) {
    fit <- stats::lm.fit(X, y)
    rss <- sum(fit$residuals^2)
    n <- length(y)
    k <- ncol(X) 
    sigma2 <- rss / n # MLE of variance (common for ICs)
    # Classical BIC/AIC: n*log(RSS/n) + k*log(n) or n*log(RSS/n) + 2k. 
    # Our n*log(sigma2) is equivalent to n*log(RSS/n).
    ic_val <- if (criterion == "BIC") {
      n * log(sigma2) + k * log(n)
    } else { 
      n * log(sigma2) + 2 * k
    }
    return(ic_val)
  }

  base_X <- matrix(1, nrow = length(y), ncol = 1)
  base_ic <- calc_ic(y, base_X)
  if (verbose) message(sprintf("[Select_SSR] Base IC (intercept-only): %.4f (criterion: %s)", base_ic, criterion))
  
  selected <- integer(0) 
  remaining <- seq_len(n_pairs) 
  current_X <- base_X 

  iter <- 0
  repeat {
    iter <- iter + 1
    if (verbose && iter > n_pairs +1 ) { message("Select_SSR: Breaking runaway loop"); break} # safety break
    if (verbose) message(sprintf("  [Select_SSR] Iteration %d, current_base_ic: %.4f, selected pairs: %s, remaining pairs: %s", 
                               iter, base_ic, paste(selected,collapse=","), paste(remaining,collapse=",")))
    best_ic_this_iter <- Inf 
    best_pair_this_iter <- NA 

    for (idx in remaining) {
      cols <- (2 * (idx - 1) + 1):(2 * idx)
      X_tmp <- cbind(current_X, U_candidates[, cols, drop = FALSE])
      ic_val <- calc_ic(y, X_tmp)
      if (verbose) message(sprintf("    Trying pair %d (cols %d-%d), IC = %.4f", idx, cols[1], cols[2], ic_val))
      if (ic_val < best_ic_this_iter) {
        best_ic_this_iter <- ic_val
        best_pair_this_iter <- idx
      }
    }

    if (is.na(best_pair_this_iter)) {
        if (verbose) message("    No best pair found in this iteration (all remaining pairs resulted in non-improvement or Inf IC).")
        break 
    }

    if (verbose) message(sprintf("    Iteration %d: Best pair to add is %d with IC = %.4f (improvement from %.4f is %.4f)", 
                               iter, best_pair_this_iter, best_ic_this_iter, base_ic, base_ic - best_ic_this_iter))

    if ((base_ic - best_ic_this_iter) > delta_threshold) { 
      cols_to_add <- (2 * (best_pair_this_iter - 1) + 1):(2 * best_pair_this_iter)
      current_X <- cbind(current_X, U_candidates[, cols_to_add, drop = FALSE])
      base_ic <- best_ic_this_iter 
      selected <- c(selected, best_pair_this_iter)
      remaining <- setdiff(remaining, best_pair_this_iter)
      if (verbose) message(sprintf("    Pair %d ADDED. New base_ic: %.4f. Selected so far: %s", 
                                 best_pair_this_iter, base_ic, paste(selected, collapse=",")))
      if (length(remaining) == 0) {
          if (verbose) message("    No more pairs remaining to test.")
          break 
      }
    } else {
      if (verbose) message(sprintf("    Best pair %d did not meet delta_threshold (%.4f <= %.4f). Stopping.", 
                                 best_pair_this_iter, base_ic - best_ic_this_iter, delta_threshold))
      break 
    }
  }

  if (length(selected) == 0) {
      if (verbose) message("[Select_SSR] No pairs selected. Returning 0-column matrix with attributes.")
      return(.empty_spec_matrix(length(y)))
  }
  cols_final <- unlist(lapply(selected, function(i) (2 * (i - 1) + 1):(2 * i)))
  res <- U_candidates[, cols_final, drop = FALSE]
  if (verbose) message(sprintf("[Select_SSR] Selected %d pairs. Final matrix dim: %s", length(selected), paste(dim(res),collapse="x")))
  
  hz <- attr(U_candidates, "freq_hz")
  rad <- attr(U_candidates, "freq_rad_s")
  if (!is.null(hz)) attr(res, "freq_hz") <- hz[selected]
  if (!is.null(rad)) attr(res, "freq_rad_s") <- rad[selected]
  attr(res, "selected_pair_indices") <- selected
  res
}

