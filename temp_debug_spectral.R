devtools::load_all("/Users/bbuchsbaum/code/fmridenoise")

# From test-spectral.R, nyquist_guard_factor test
set.seed(789)
TR_val <- 1.0
n_timepoints <- 200
time_vector_sec <- (seq_len(n_timepoints) - 1) * TR_val
nyquist <- 1 / (2 * TR_val) # 0.5 Hz

# New strategy: Very strong, low-frequency signal
freq_strong <- 0.1 # Hz (well below Nyquist)
signal_component <- sin(2 * pi * freq_strong * time_vector_sec) * 100 # Very high amplitude
mean_resid_data <- signal_component + rnorm(n_timepoints, sd = 0.001)   # Extremely low noise

message("--- Testing ndx_spectral_sines with nyquist_guard_factor = 1.0 (MUCH Stronger, Lower Freq Signal) ---")
spectral_regressors_unguarded <- ndx_spectral_sines(
  mean_resid_data, 
  TR = TR_val, 
  n_sine_candidates = 1, 
  nyquist_guard_factor = 1.0,
  selection_criterion = "AIC",
  selection_delta_threshold = 0.1,
  verbose = TRUE
)

if (is.null(spectral_regressors_unguarded)) {
  message("spectral_regressors_unguarded IS NULL")
} else {
  message(paste("spectral_regressors_unguarded dims:", paste(dim(spectral_regressors_unguarded), collapse="x")))
  message("Identified freqs (Hz):")
  print(attr(spectral_regressors_unguarded, "freq_hz"))
} 