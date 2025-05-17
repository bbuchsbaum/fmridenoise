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
  signal_near_nyquist <- sin(2 * pi * freq_near_nyquist * time_vector_sec)
  mean_resid_data <- signal_near_nyquist + rnorm(n_timepoints, sd = 0.1)

  # With default guard factor (0.9), freq_near_nyquist (0.475Hz) should be excluded 
  # as 0.9 * 0.5 = 0.45 Hz is the upper limit.
  spectral_regressors_guarded <- ndx_spectral_sines(mean_resid_data, TR = TR_val, n_sine_candidates = 1, nyquist_guard_factor = 0.9)
  if (!is.null(spectral_regressors_guarded)) {
    identified_freqs_guarded <- attr(spectral_regressors_guarded, "freq_hz")
    expect_true(all(identified_freqs_guarded < freq_near_nyquist), 
                "Identified frequencies should be less than the near-Nyquist signal if guarded.")
  } else {
    expect_null(spectral_regressors_guarded, "Expected no peaks or NULL due to guard factor excluding the main peak.")
  }

  # With guard factor = 1.0 (or > 0.95), it should be found.
  spectral_regressors_unguarded <- ndx_spectral_sines(mean_resid_data, TR = TR_val, n_sine_candidates = 1, nyquist_guard_factor = 1.0)
  expect_true(!is.null(spectral_regressors_unguarded), "Should find the peak when guard factor is 1.0.")
  if (!is.null(spectral_regressors_unguarded)) {
    identified_freqs_unguarded <- attr(spectral_regressors_unguarded, "freq_hz")
    expect_true(any(abs(identified_freqs_unguarded - freq_near_nyquist) < 0.01), 
                "The near-Nyquist frequency should be identified when guard factor is 1.0.")
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
    "k_tapers became 0 after adjustment"
  )
  expect_null(res)
}) 