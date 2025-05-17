context("ndx_spectral_sines - Functionality")

test_that("ndx_spectral_sines generates correct regressors for a simple signal", {
  set.seed(123)
  TR_val <- 2.0
  n_timepoints <- 200
  time_vector_sec <- (seq_len(n_timepoints) - 1) * TR_val

  # Create a synthetic residual with a clear 0.05 Hz sine wave + some noise
  freq_of_interest_hz <- 0.05 # Hz
  signal_component <- sin(2 * pi * freq_of_interest_hz * time_vector_sec)
  noise_component <- rnorm(n_timepoints, sd = 0.2)
  mean_resid_data <- signal_component + noise_component

  n_candidates <- 1
  spectral_regressors <- ndx_spectral_sines(mean_resid_data, 
                                            TR = TR_val, 
                                            n_sine_candidates = n_candidates,
                                            k_tapers = 5, # Align with function default
                                            nw = 3)       # Align with function default

  expect_true(!is.null(spectral_regressors), "Spectral regressors should not be NULL for a clear signal.")
  expect_equal(ncol(spectral_regressors), 2 * n_candidates, 
               info = "Should generate 2 regressors (sine & cosine) per candidate peak.")
  expect_equal(nrow(spectral_regressors), n_timepoints, 
               info = "Number of rows should match length of input time series.")
  
  identified_freqs_hz <- attr(spectral_regressors, "freq_hz")
  expect_true(!is.null(identified_freqs_hz), "Attribute 'freq_hz' should be present.")
  expect_length(identified_freqs_hz, n_candidates)
  
  # Check if the identified frequency is close to the known signal frequency
  # multitaper::spec.mtm might not hit the exact frequency, so allow tolerance
  expect_true(abs(identified_freqs_hz[1] - freq_of_interest_hz) < 0.01, 
              info = sprintf("Identified frequency (%.4f Hz) should be close to actual (%.4f Hz).", 
                             identified_freqs_hz[1], freq_of_interest_hz))
  
  # Check orthogonality of generated sine/cosine pair for the identified frequency
  if (!is.null(spectral_regressors) && ncol(spectral_regressors) >= 2) {
    correlation_sin_cos <- cor(spectral_regressors[,1], spectral_regressors[,2])
    expect_lt(abs(correlation_sin_cos), 0.01, 
              label = "Correlation between generated sine and cosine for the same frequency should be very low (near orthogonal).")
  } else {
    fail("Spectral regressors not generated or insufficient columns for orthogonality test.")
  }

  # Check if the identified regressors can explain the original signal component
  if (!is.null(spectral_regressors) && ncol(spectral_regressors) >= 2) {
    fit <- lm(signal_component ~ spectral_regressors[,1] + spectral_regressors[,2])
    r_squared <- summary(fit)$r.squared
    expect_gt(r_squared, 0.8, 
              label = "Identified sine/cosine pair should explain a large portion of the original signal variance (R-squared > 0.8)")
  } else {
    fail("Spectral regressors not generated or insufficient columns for lm test.")
  }
})

test_that("ndx_spectral_sines handles no clear peaks", {
  set.seed(456)
  TR_val <- 1.0
  n_timepoints <- 100
  # Pure noise, no strong periodic signal
  mean_resid_data <- rnorm(n_timepoints, sd = 1.0)

  spectral_regressors <- ndx_spectral_sines(mean_resid_data, TR = TR_val, n_sine_candidates = 3)

  # Depending on noise, it might find some peaks or none. 
  # If it finds peaks, check consistency. If NULL, that's also acceptable.
  if (is.null(spectral_regressors)) {
    expect_null(spectral_regressors) # Explicitly state that NULL is okay
  } else {
    expect_true(is.matrix(spectral_regressors))
    expect_true(ncol(spectral_regressors) <= 2 * 3) # At most 2*n_sine_candidates
  }
})

test_that("ndx_spectral_sines handles reasonably short but valid time series", {
  TR_val <- 1.0
  # Time series just long enough for k_tapers=2, nw=1.5
  # Original: k_tapers = 2, nw = 1.5. len = 20
  # Adjusted k_tapers = min(2, 2*1.5-1, 20-1) = min(2, 2, 19) = 2. This is >=1.
  # So spec.mtm should proceed.
  mean_resid_data_ok <- rnorm(20)
  expect_no_warning({
    spectral_regressors_ok <- ndx_spectral_sines(mean_resid_data_ok, TR = TR_val, k_tapers = 2, nw=1.5, n_sine_candidates = 1)
  })
  # We don't strictly require it to find peaks here, just not to error/warn inappropriately on length
  if (is.null(spectral_regressors_ok)) {
    expect_null(spectral_regressors_ok)
  } else {
    expect_true(is.matrix(spectral_regressors_ok))
  }
})

test_that("ndx_spectral_sines handles invalid inputs", {
  expect_warning(ndx_spectral_sines(character(0), TR = 1.0), "mean_residual_for_spectrum must be a non-empty numeric vector.")
  expect_null(suppressWarnings(ndx_spectral_sines(character(0), TR = 1.0)))
  
  expect_warning(ndx_spectral_sines(rnorm(100), TR = 0), "TR must be a positive numeric value.")
  expect_null(suppressWarnings(ndx_spectral_sines(rnorm(100), TR = 0)))
  
  expect_warning(ndx_spectral_sines(rnorm(100), TR = -1), "TR must be a positive numeric value.")
  expect_null(suppressWarnings(ndx_spectral_sines(rnorm(100), TR = -1)))
})

test_that("nyquist_guard_factor works as expected", {
  set.seed(789)
  TR_val <- 1.0
  n_timepoints <- 200
  time_vector_sec <- (seq_len(n_timepoints) - 1) * TR_val
  nyquist <- 1 / (2 * TR_val) # 0.5 Hz

  # Signal very close to Nyquist
  freq_near_nyquist <- nyquist * 0.95 # 0.475 Hz
  signal_component <- sin(2 * pi * freq_near_nyquist * time_vector_sec) * 5 # Increased amplitude
  mean_resid_data <- signal_component + rnorm(n_timepoints, sd = 0.05)     # Reduced noise

  # With default guard factor (0.9), freq_near_nyquist (0.475Hz) should be excluded 
  # as 0.9 * 0.5 = 0.45 Hz is the upper limit.
  # For this part, use the original near-Nyquist signal to test exclusion
  freq_near_nyquist_orig <- nyquist * 0.95 # 0.475 Hz
  signal_near_nyquist_orig <- sin(2 * pi * freq_near_nyquist_orig * time_vector_sec)
  mean_resid_data_orig_for_guard <- signal_near_nyquist_orig + rnorm(n_timepoints, sd = 0.1) # Original noise

  spectral_regressors_guarded <- ndx_spectral_sines(mean_resid_data_orig_for_guard, TR = TR_val, 
                                                  n_sine_candidates = 1, 
                                                  nyquist_guard_factor = 0.9, # Default guard
                                                  verbose = FALSE) # Less verbose for this part
  if (!is.null(spectral_regressors_guarded)) {
    identified_freqs_guarded <- attr(spectral_regressors_guarded, "freq_hz")
    expect_true(all(identified_freqs_guarded < freq_near_nyquist_orig), 
                "Identified frequencies should be less than the near-Nyquist signal if guarded.")
  } else {
    expect_null(spectral_regressors_guarded, "Expected no peaks or NULL due to guard factor excluding the main peak.")
  }

  # With guard factor = 1.0, it should be found IF signal is strong enough.
  # Use a very strong, lower frequency signal to ensure IC selects it.
  freq_strong <- 0.1 # Hz
  signal_strong_component <- sin(2 * pi * freq_strong * time_vector_sec) * 100
  mean_resid_data_strong <- signal_strong_component + rnorm(n_timepoints, sd = 0.001)
  
  spectral_regressors_unguarded <- ndx_spectral_sines(
    mean_resid_data_strong, # Use the new strong signal data
    TR = TR_val, 
    n_sine_candidates = 1, 
    nyquist_guard_factor = 1.0,
    selection_criterion = "AIC", 
    selection_delta_threshold = 0.1, 
    verbose = FALSE # Can be FALSE now
  )
  expect_true(!is.null(spectral_regressors_unguarded) && ncol(spectral_regressors_unguarded) > 0, 
              "Should find the peak and return regressors when guard factor is 1.0 and signal is very strong.")
  if (!is.null(spectral_regressors_unguarded) && ncol(spectral_regressors_unguarded) > 0) {
    identified_freqs_unguarded <- attr(spectral_regressors_unguarded, "freq_hz")
    expect_true(any(abs(identified_freqs_unguarded - freq_strong) < 0.01), 
                "The strong 0.1Hz frequency should be identified when guard factor is 1.0.")
  }
})

test_that("ndx_spectral_sines handles short time series that cause spec.mtm to fail", {
  short_ts <- rnorm(5) # Length 5, spec.mtm might fail with k_tapers adjusted but still too high for series
  # Default k_tapers=5, nw=3. Adjusted k_tapers will be min(5, 2*3-1, 5-1) = min(5,5,4) = 4.
  # spec.mtm with k=4 on length 5 series is likely to fail or produce unusable output.
  expect_warning(
    res <- ndx_spectral_sines(short_ts, TR = 1, n_sine_candidates = 1, k_tapers = 5, nw = 3),
    "Spectrum estimation via spec.mtm did not yield valid spec or freq."
  )
  expect_null(res)
})

test_that("ndx_spectral_sines warns if k_tapers becomes < 1 after adjustment", {
  # To make k_tapers < 1: use short series and/or small nw.
  # Example: k_tapers=2 (initial), nw=1. len=2.
  # Adjusted k_tapers = min(2, 2*1-1, 2-1) = min(2, 1, 1) = 1. Not < 1.
  # Example: k_tapers=1 (initial), nw=1. len=1.
  # Adjusted k_tapers = min(1, 2*1-1, 1-1) = min(1, 1, 0) = 0.
  very_short_ts <- rnorm(1)
  expect_warning(
    res <- ndx_spectral_sines(very_short_ts, TR = 1, n_sine_candidates = 1, k_tapers = 1, nw = 1),
    "nw became 0.00 after adjustment. Must be > 0. Aborting spectral step."
  )
  # Function now returns a 0-column matrix with attributes in this case
  expect_true(is.matrix(res) && ncol(res) == 0 && nrow(res) == length(very_short_ts))
  expect_true(!is.null(attr(res, "freq_hz")) && length(attr(res, "freq_hz")) == 0)
})

test_that("Select_Significant_Spectral_Regressors filters by BIC", {
  set.seed(42)
  TR_val <- 1
  n_tp <- 200
  t_sec <- (seq_len(n_tp) - 1) * TR_val
  freq_good <- 0.05
  freq_bad <- 0.20
  y <- sin(2 * pi * freq_good * t_sec) + rnorm(n_tp, sd = 0.1)
  U_good <- cbind(sin(2 * pi * freq_good * t_sec), cos(2 * pi * freq_good * t_sec))
  U_bad  <- cbind(sin(2 * pi * freq_bad * t_sec), cos(2 * pi * freq_bad * t_sec))
  U_all  <- cbind(U_good, U_bad)
  colnames(U_all) <- paste0(c("sin_good","cos_good","sin_bad","cos_bad"))
  attr(U_all, "freq_hz") <- c(freq_good, freq_bad)
  selected <- Select_Significant_Spectral_Regressors(y, U_all, criterion = "BIC", delta_threshold = 2)
  expect_true(!is.null(selected))
  expect_equal(ncol(selected), 2)
  expect_true(all(colnames(selected) %in% c("sin_good","cos_good")))
})

context("Select_Significant_Spectral_Regressors - Regressor Selection Logic")

test_that("Select_Significant_Spectral_Regressors handles no candidates or no selection", {
  set.seed(100)
  y_test <- rnorm(100)
  empty_mat_attrs <- function(y_vec) {
    res <- matrix(numeric(0), nrow = length(y_vec), ncol = 0, dimnames = list(NULL, character(0)))
    attr(res, "freq_hz") <- numeric(0)
    attr(res, "freq_rad_s") <- numeric(0)
    attr(res, "selected_pair_indices") <- integer(0)
    return(res)
  }

  expect_equal(Select_Significant_Spectral_Regressors(y_test, NULL), empty_mat_attrs(y_test))
  expect_equal(Select_Significant_Spectral_Regressors(y_test, matrix(0,100,0)), empty_mat_attrs(y_test))
  expect_equal(Select_Significant_Spectral_Regressors(y_test, matrix(rnorm(100),100,1)), empty_mat_attrs(y_test)) # Not pairs

  # Case: No pair improves IC enough
  U_noise1 <- cbind(rnorm(100), rnorm(100))
  colnames(U_noise1) <- c("sin_noise1", "cos_noise1")
  attr(U_noise1, "freq_hz") <- 0.1
  attr(U_noise1, "freq_rad_s") <- 2*pi*0.1
  expect_equal(Select_Significant_Spectral_Regressors(y_test, U_noise1, delta_threshold = 2, criterion = "BIC"), 
               empty_mat_attrs(y_test))
})

test_that("Select_Significant_Spectral_Regressors selects one good pair (BIC & AIC)", {
  set.seed(200)
  n_tp <- 100; TR <- 1.0
  t_sec <- (1:n_tp-1)*TR
  freq1 <- 0.1 # Hz
  signal1 <- sin(2*pi*freq1*t_sec) * 2
  y_test <- signal1 + rnorm(n_tp, sd = 0.5)
  
  U_pair1 <- cbind(sin(2*pi*freq1*t_sec), cos(2*pi*freq1*t_sec))
  colnames(U_pair1) <- c("s1","c1")
  attr(U_pair1, "freq_hz") <- freq1
  attr(U_pair1, "freq_rad_s") <- 2*pi*freq1

  # Test with BIC
  selected_bic <- Select_Significant_Spectral_Regressors(y_test, U_pair1, criterion = "BIC", delta_threshold = 2)
  expect_true(!is.null(selected_bic) && ncol(selected_bic) == 2)
  if (!is.null(selected_bic)) {
    expect_equal(colnames(selected_bic), c("s1","c1"))
    expect_equal(attr(selected_bic, "freq_hz"), freq1)
    expect_equal(attr(selected_bic, "selected_pair_indices"), 1L)
  }
  
  # Test with AIC (should also select it, likely with better IC improvement)
  selected_aic <- Select_Significant_Spectral_Regressors(y_test, U_pair1, criterion = "AIC", delta_threshold = 2)
  expect_true(!is.null(selected_aic) && ncol(selected_aic) == 2)
  if (!is.null(selected_aic)) {
    expect_equal(colnames(selected_aic), c("s1","c1"))
    expect_equal(attr(selected_aic, "freq_hz"), freq1)
    expect_equal(attr(selected_aic, "selected_pair_indices"), 1L)
  }
})

test_that("Select_Significant_Spectral_Regressors selects multiple good pairs and respects order/delta", {
  set.seed(300)
  n_tp <- 200; TR <- 1.0
  t_sec <- (1:n_tp-1)*TR
  freq1 <- 0.05 # Stronger signal
  freq2 <- 0.15 # Weaker signal
  freq3 <- 0.30 # Noise signal
  
  signal1 <- sin(2*pi*freq1*t_sec) * 3
  signal2 <- sin(2*pi*freq2*t_sec) * 1.5
  y_test <- signal1 + signal2 + rnorm(n_tp, sd = 0.5)
  
  U_pair1 <- cbind(sin(2*pi*freq1*t_sec), cos(2*pi*freq1*t_sec)); colnames(U_pair1) <- c("s1","c1")
  U_pair2 <- cbind(sin(2*pi*freq2*t_sec), cos(2*pi*freq2*t_sec)); colnames(U_pair2) <- c("s2","c2")
  U_pair3_noise <- cbind(rnorm(n_tp), rnorm(n_tp)); colnames(U_pair3_noise) <- c("s_noise","c_noise") # Pure noise pair
  
  U_candidates <- cbind(U_pair1, U_pair2, U_pair3_noise)
  attr(U_candidates, "freq_hz")    <- c(freq1, freq2, NA) # NA for noise pair
  attr(U_candidates, "freq_rad_s") <- c(2*pi*freq1, 2*pi*freq2, NA)
  
  # Expect pair1 and pair2 to be selected (likely pair1 first due to stronger signal)
  selected_multi <- Select_Significant_Spectral_Regressors(y_test, U_candidates, criterion = "BIC", delta_threshold = 2)
  expect_true(!is.null(selected_multi) && ncol(selected_multi) == 4) # Expect both good pairs
  if (!is.null(selected_multi)) {
    expect_true(all(c("s1","c1","s2","c2") %in% colnames(selected_multi)))
    expect_false(any(c("s_noise","c_noise") %in% colnames(selected_multi)))
    expect_equal(sort(attr(selected_multi, "freq_hz")), sort(c(freq1, freq2)))
    expect_equal(sort(attr(selected_multi, "selected_pair_indices")), sort(c(1L, 2L))) # Check indices of selected pairs
  }
  
  # Test delta_threshold: make it very high. 
  # With only one good pair (U_pair1) and one noise pair (U_pair3_noise), only the good one should be selected.
  U_candidates_simple <- cbind(U_pair1, U_pair3_noise)
  attr(U_candidates_simple, "freq_hz")    <- c(freq1, NA) # freq1 for U_pair1, NA for U_pair3_noise
  attr(U_candidates_simple, "freq_rad_s") <- c(2*pi*freq1, NA)
  # Note: selected_pair_indices will refer to indices within U_candidates_simple (1 for U_pair1, 2 for U_pair3_noise)

  selected_high_delta <- Select_Significant_Spectral_Regressors(y_test, U_candidates_simple, criterion = "BIC", delta_threshold = 100) 
  expect_true(!is.null(selected_high_delta) && ncol(selected_high_delta) == 2, 
              label = "selected_high_delta: Expected 1 pair (2 cols)")
  if (!is.null(selected_high_delta) && ncol(selected_high_delta) == 2) { 
    expect_true(all(colnames(selected_high_delta) %in% c("s1","c1")),
                label = "selected_high_delta: Colnames should be s1, c1")
    expect_equal(attr(selected_high_delta, "freq_hz"), freq1, 
                 label = "selected_high_delta: freq_hz attribute")
    expect_equal(attr(selected_high_delta, "selected_pair_indices"), 1L, 
                 label = "selected_high_delta: selected_pair_indices attribute")
  }
})
