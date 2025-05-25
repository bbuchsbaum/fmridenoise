context("ndx_spectral_sines - environment cleanup")

test_that("mt_res is reset locally without affecting global env", {
  old_exists <- exists("mt_res", envir = globalenv(), inherits = FALSE)
  if (old_exists) old_val <- get("mt_res", envir = globalenv())

  assign("mt_res", "global_value", envir = globalenv())
  short_ts <- rnorm(5)

  expect_warning(
    res <- ndx_spectral_sines(short_ts, TR = 1, n_sine_candidates = 1, k_tapers = 5, nw = 3),
    "Spectrum estimation via spec.mtm did not yield valid spec or freq."
  )
  expect_null(res)
  expect_identical(get("mt_res", envir = globalenv()), "global_value")

  if (old_exists) {
    assign("mt_res", old_val, envir = globalenv())
  } else {
    rm("mt_res", envir = globalenv())
  }
})
