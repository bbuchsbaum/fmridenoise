context("diagnostic plot helpers")

# sample data for plots
sample_diagnostics <- list(
  list(DES = 0.1, ljung_box_p = 0.3),
  list(DES = 0.2, ljung_box_p = 0.6)
)

sample_betas <- list(
  list(matrix(rnorm(4), 2, 2), matrix(rnorm(4), 2, 2)),
  list(matrix(rnorm(4), 2, 2), matrix(rnorm(4), 2, 2))
)

sample_pass0 <- matrix(rnorm(20), 10, 2)
sample_final <- matrix(rnorm(20), 10, 2)

sample_S <- matrix(rnorm(20), 10, 2)

sample_gdlite <- list(r2_vals_by_k = c(0.1, 0.3, 0.4), optimal_k = 2)

TR_test <- 2

# ---- Successful cases ----

test_that("plot helpers create pngs and return paths", {
  outdir <- tempfile("plots")
  dir.create(outdir)

  des_path <- ndx:::ndx_plot_des_per_pass(sample_diagnostics, outdir)
  expect_true(file.exists(des_path))

  beta_path <- ndx:::ndx_plot_beta_stability(sample_betas, outdir)
  expect_true(file.exists(beta_path))

  psd_path <- ndx:::ndx_plot_residual_psd(sample_pass0, sample_final, TR_test, outdir)
  expect_true(file.exists(psd_path))

  ljung_path <- ndx:::ndx_plot_ljung_box_pvalues(sample_diagnostics, outdir)
  expect_true(file.exists(ljung_path))

  carpet_path <- ndx:::ndx_plot_spike_carpet(sample_S, outdir)
  expect_true(file.exists(carpet_path))

  gdlite_path <- ndx:::ndx_plot_gdlite_r2_by_k(sample_gdlite, outdir)
  expect_true(file.exists(gdlite_path))
})

# ---- Edge cases ----

test_that("plot helpers handle invalid input", {
  expect_error(ndx:::ndx_plot_des_per_pass(list(), tempdir()))
  expect_error(ndx:::ndx_plot_beta_stability(list(), tempdir()))
  expect_error(ndx:::ndx_plot_residual_psd(matrix(1,2,2), matrix(1,3,2), TR_test, tempdir()))
  expect_null(ndx:::ndx_plot_ljung_box_pvalues(list(list(ljung_box_p = NA_real_)), tempdir()))
  expect_error(ndx:::ndx_plot_spike_carpet(1:5, tempdir()))
  expect_error(ndx:::ndx_plot_gdlite_r2_by_k(list(), tempdir()))
})
