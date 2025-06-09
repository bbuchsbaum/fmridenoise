context("progressive viz internal helpers")

# Test ndx_calculate_annihilation_verdict for a high variance ratio

test_that("ndx_calculate_annihilation_verdict produces expected html", {
  wf <- list(
    gdlite_pcs = matrix(1, 1, 1),
    rpca_orthogonalized = matrix(sqrt(1.5), 1, 1),
    spectral_orthogonalized = NULL
  )
  html <- ndx:::ndx_calculate_annihilation_verdict(wf)
  expect_true(grepl("Annihilation", html, fixed = TRUE))
  expect_true(grepl("#d9534f", html, fixed = TRUE))
})

# When required components are missing, should return empty string

test_that("ndx_calculate_annihilation_verdict handles missing components", {
  expect_identical(ndx:::ndx_calculate_annihilation_verdict(list()), "")
})

# Test create_des_progression_plot with valid input and error case

test_that("create_des_progression_plot returns plotly object", {
  wf <- list(diagnostics_per_pass = list(list(DES = 0.1), list(DES = 0.2)))
  p <- ndx:::create_des_progression_plot(wf, c("S1", "S2"))
  expect_s3_class(p, "plotly")
})

test_that("create_des_progression_plot errors with missing diagnostics", {
  expect_error(ndx:::create_des_progression_plot(list(), c("S1")))
})

# Test create_psd_progression_plot

test_that("create_psd_progression_plot returns plotly object", {
  wf <- list(
    pass0_residuals = matrix(rnorm(20), nrow = 20, ncol = 1),
    Y_residuals_final_unwhitened = matrix(rnorm(20), nrow = 20, ncol = 1)
  )
  p <- ndx:::create_psd_progression_plot(wf, c("S1", "S4"))
  expect_s3_class(p, "plotly")
})

