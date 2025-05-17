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
#'
#' @return A matrix with `2 * n_selected_peaks` columns and `length(mean_residual_for_spectrum)`
#'   rows. Each pair of columns represents sine and cosine regressors for an
#'   identified peak frequency. The attribute "freq" contains the frequencies
#'   (in Hz) of the selected peaks. Returns NULL if no peaks are found or if
#'   input is unsuitable.
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
#'     print(attr(U_sines, "freq"))
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
                               k_tapers = 5, nw = 3) {

  if (!is.numeric(mean_residual_for_spectrum) || length(mean_residual_for_spectrum) == 0) {
    warning("mean_residual_for_spectrum must be a non-empty numeric vector.")
    return(NULL)
  }
  if (!is.numeric(TR) || TR <= 0) {
    warning("TR must be a positive numeric value.")
    return(NULL)
  }

  # Adjust k_tapers based on nw and series length
  # Ensure k_tapers is at least 1 and respects multitaper constraints.
  # length(mean_residual_for_spectrum) - 1 is to ensure k_tapers is strictly less than N for some spec.mtm internals or typical usage.
  # 2*nw - 1 is a constraint from the multitaper theory.
  k_tapers <- min(k_tapers, 2 * nw - 1, length(mean_residual_for_spectrum) - 1)
  if (k_tapers < 1) {
    warning(sprintf("k_tapers became %d after adjustment (nw=%s, series_length=%d). Must be >= 1. Aborting spectral step.", 
                    as.integer(k_tapers), nw, length(mean_residual_for_spectrum)))
    return(NULL)
  }

  mt_res <- NULL
  tryCatch({
    mt_res <- multitaper::spec.mtm(mean_residual_for_spectrum,
                                   k       = k_tapers,
                                   nw      = nw,
                                   deltat  = TR,
                                   jackknife = FALSE,
                                   returnInternals = FALSE,
                                   plot    = FALSE)
  }, error = function(e) {
    warning(paste("multitaper::spec.mtm failed:", e$message))
    mt_res <<- NULL
  })

  if (is.null(mt_res) || is.null(mt_res$spec) || is.null(mt_res$freq) || length(mt_res$spec) == 0) {
    warning("Spectrum estimation via spec.mtm did not yield valid spec or freq.")
    return(NULL)
  }

  nyquist_freq  <- 1 / (2 * TR)
  # Ensure guard factor is within reasonable bounds
  guard_factor <- max(0, min(1, nyquist_guard_factor)) 
  
  freq_upper_bound <- guard_factor * nyquist_freq
  keep_indices <- mt_res$freq > 0 & mt_res$freq <= freq_upper_bound # Exclude DC, respect guard

  if (sum(keep_indices) == 0) {
    warning("No frequencies to search for peaks after applying guard factor and excluding DC.")
    return(NULL)
  }

  spec_to_search <- mt_res$spec[keep_indices]
  freq_to_search <- mt_res$freq[keep_indices]

  # findpeaks returns a matrix: col1=amplitude, col2=peak_index, col3=start_index, col4=end_index
  # We need indices relative to spec_to_search
  peak_info_all <- pracma::findpeaks(spec_to_search, nups = 1, ndowns = 1, sortstr = TRUE)

  if (is.null(peak_info_all) || nrow(peak_info_all) == 0) {
    message("No initial peaks found in the spectrum by pracma::findpeaks.")
    return(NULL)
  }
  
  # Peak prominence filter
  prom_thresh <- median(spec_to_search, na.rm = TRUE) + 3 * mad(spec_to_search, na.rm = TRUE)
  # Ensure threshold is a single finite number, if mad is 0 or spec_to_search is flat, this could be tricky
  if (!is.finite(prom_thresh)) {
    prom_thresh <- -Inf # Effectively disable filter if MAD is zero or issues
    warning("Prominence threshold for peak picking was not finite (e.g. MAD was zero). Filter effectively disabled.")
  }
  
  # Keep peaks with amplitude (col 1) > prom_thresh
  # peak_info_all is already sorted by amplitude (col 1) due to sortstr=TRUE
  peak_info <- peak_info_all[peak_info_all[,1] > prom_thresh, , drop = FALSE]

  if (is.null(peak_info) || nrow(peak_info) == 0) {
    message(sprintf("No significant peaks found after prominence filter (threshold: %.4g).", prom_thresh))
    return(NULL)
  }
  
  # Select top n_sine_candidates peaks based on their amplitude (peak_info already sorted by findpeaks with sortstr=TRUE)
  num_peaks_to_select <- min(n_sine_candidates, nrow(peak_info))
  selected_peak_indices_in_searched_spec <- peak_info[1:num_peaks_to_select, 2]

  # Frequencies corresponding to these selected peaks
  selected_frequencies_hz <- freq_to_search[selected_peak_indices_in_searched_spec]
  
  if (length(selected_frequencies_hz) == 0) {
    message("No peaks selected after filtering.")
    return(NULL)
  }

  omega_rad_s <- 2 * pi * selected_frequencies_hz # Convert Hz to radians/sec
  time_vector_sec <- (seq_along(mean_residual_for_spectrum) - 1) * TR # Time vector in seconds, starting at 0

  U_sin_raw <- sapply(omega_rad_s, function(w) sin(w * time_vector_sec))
  U_cos_raw <- sapply(omega_rad_s, function(w) cos(w * time_vector_sec))
  
  # Ensure matrix even if only one frequency
  if (is.vector(U_sin_raw)) U_sin_raw <- matrix(U_sin_raw, ncol = 1)
  if (is.vector(U_cos_raw)) U_cos_raw <- matrix(U_cos_raw, ncol = 1)
  
  # Normalize columns to have unit norm (for exact orthogonality)
  # center = FALSE because we don't want to de-mean sines/cosines
  # scale = sqrt(colSums(X^2)) to divide by L2 norm
  U_sin <- apply(U_sin_raw, 2, function(col) {
    norm_val <- sqrt(sum(col^2))
    if (norm_val > .Machine$double.eps) col / norm_val else col
  })
  U_cos <- apply(U_cos_raw, 2, function(col) {
    norm_val <- sqrt(sum(col^2))
    if (norm_val > .Machine$double.eps) col / norm_val else col
  })
  
  # Ensure matrix structure is preserved after apply if only one column
  if (is.vector(U_sin)) U_sin <- matrix(U_sin, ncol = 1)
  if (is.vector(U_cos)) U_cos <- matrix(U_cos, ncol = 1)

  # Ensure column names are unique if multiple frequencies are identical (should be rare with findpeaks)
  colnames(U_sin) <- paste0("sin_f", sprintf("%.4f", selected_frequencies_hz))
  colnames(U_cos) <- paste0("cos_f", sprintf("%.4f", selected_frequencies_hz))
  
  U_spectral_sines <- cbind(U_sin, U_cos)

  attr(U_spectral_sines, "freq_hz") <- selected_frequencies_hz
  attr(U_spectral_sines, "freq_rad_s") <- omega_rad_s
  
  message(sprintf("Generated %d sine/cosine pairs from %d spectral peaks.", 
                  ncol(U_spectral_sines) / 2, length(selected_frequencies_hz)))
  return(U_spectral_sines)
} 