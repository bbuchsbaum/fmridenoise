context("ndx_generate_html_report")

TR_test <- 2.0
n_time <- 50
n_vox <- 4
Y0 <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
pass0_res <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
workflow_mock <- list(
  Y_residuals_final_unwhitened = Y0,
  diagnostics_per_pass = list(
    list(DES = 0.1),
    list(
      DES = 0.2,
      k_rpca_global = 3,
      num_hrf_clusters = 1,
      num_spectral_sines = 1,
      lambda_parallel_noise = 0.5,
      lambda_perp_signal = 0.1,
      rho_noise_projection = 0.25,
      ljung_box_p = 0.6
    )
  ),
  betas_per_pass = list(
    matrix(rnorm(n_vox * n_vox), n_vox, n_vox)
  ),
  beta_history_per_pass = list(
    matrix(rnorm(n_time * n_vox), n_time, n_vox),
    matrix(rnorm(n_time * n_vox), n_time, n_vox)
  ),
  num_passes_completed = 2,
  S_matrix_rpca_final = matrix(0, n_time, n_vox)
)

test_that("HTML report is generated", {
  tmpdir <- tempfile("diag")
  dir.create(tmpdir)
  expect_no_error({
    html_path <- ndx_generate_html_report(workflow_mock, pass0_res, TR_test, output_dir = tmpdir)
  })
  expect_true(file.exists(file.path(tmpdir, "ndx_diagnostic_report.html")))
  expect_true(file.exists(file.path(tmpdir, "beta_stability_per_pass.png")))
  expect_true(file.exists(file.path(tmpdir, "ljung_box_pvalues.png")))
})

test_that("JSON certificate is generated with expected fields", {
  tmpjson <- tempfile("cert", fileext = ".json")
  expect_no_error({
    json_path <- ndx_generate_json_certificate(workflow_mock, output_path = tmpjson)
  })
  expect_true(file.exists(tmpjson))
  cert <- jsonlite::read_json(tmpjson)
  expect_true(all(c(
    "ndx_version",
    "final_DES",
    "passes_converged",
    "ljung_box_p",
    "var_ratio",
    "verdict",
    "final_rho_noise_projection",
    "final_adaptive_hyperparameters"
  ) %in% names(cert)))
  expect_true(all(c(
    "k_rpca_global",
    "num_hrf_clusters",
    "num_spectral_sines",
    "lambda_parallel_noise",
    "lambda_perp_signal"
  ) %in% names(cert$final_adaptive_hyperparameters)))
})
