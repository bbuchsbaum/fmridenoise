context("ndx_spectral_sines helpers")

test_that(".validate_spectral_inputs centers data and enforces defaults", {
  vec <- 1:10
  res <- fmridenoise:::`.validate_spectral_inputs`(vec, TR = 1,
                                                   n_sine_candidates = 2,
                                                   selection_delta_threshold = 3)
  expect_null(res$error)
  expect_equal(sum(res$vec), 0)
  expect_equal(res$n_sine_candidates, 2)
  expect_equal(res$selection_delta_threshold, 3)

  expect_warning(
    res_bad <- fmridenoise:::`.validate_spectral_inputs`(vec, TR = 1,
                                                        n_sine_candidates = -1,
                                                        selection_delta_threshold = 1),
    "n_sine_candidates"
  )
  expect_equal(res_bad$n_sine_candidates, 6)
})

test_that(".adjust_mtm_params constrains k_tapers and nw", {
  adj <- fmridenoise:::`.adjust_mtm_params`(k_tapers = 10, nw = 2,
                                            series_length = 20, verbose = FALSE)
  expect_null(adj$error)
  expect_equal(adj$nw, 2)
  expect_equal(adj$k_tapers, 3)

  adj_bad <- fmridenoise:::`.adjust_mtm_params`(k_tapers = 2, nw = 1,
                                                series_length = 1, verbose = FALSE)
  expect_true(is.matrix(adj_bad$error))
  expect_equal(nrow(adj_bad$error), 1)
})

test_that(".estimate_spectrum produces spectrum or NULL on failure", {
  ts <- sin(2 * pi * 0.1 * (0:99))
  est <- fmridenoise:::`.estimate_spectrum`(ts, TR = 1, k_tapers = 3, nw = 2)
  expect_true(!is.null(est))
  expect_true(all(c("spec", "freq") %in% names(est)))

  expect_warning(
    est_bad <- fmridenoise:::`.estimate_spectrum`(rnorm(5), TR = 1,
                                                 k_tapers = 5, nw = 3),
    "spec.mtm"
  )
  expect_null(est_bad)
})

test_that(".find_spectral_peaks detects correct frequency", {
  ts <- sin(2 * pi * 0.1 * (0:199))
  est <- fmridenoise:::`.estimate_spectrum`(ts, TR = 1, k_tapers = 5, nw = 3)
  freqs <- fmridenoise:::`.find_spectral_peaks`(est, TR = 1,
                                                nyquist_guard_factor = 1,
                                                n_sine_candidates = 1,
                                                verbose = FALSE)
  expect_true(length(freqs) >= 1)
  expect_lt(abs(freqs[1] - 0.1), 0.01)
})

test_that(".build_sine_regressors constructs proper matrix", {
  U <- fmridenoise:::`.build_sine_regressors`(0.1, series_length = 10, TR = 1)
  expect_equal(dim(U), c(10, 2))
  expect_equal(attr(U, "freq_hz"), 0.1)
  expect_equal(colnames(U), c("sin_f0.1000", "cos_f0.1000"))
})
