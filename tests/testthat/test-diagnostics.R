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
      num_spectral_sines = 1,
      lambda_parallel_noise = 0.5,
      lambda_perp_signal = 0.1,
      rho_noise_projection = 0.25
    )
  ),
  beta_history_per_pass = list(
    list(matrix(rnorm(n_vox), 1, n_vox), matrix(rnorm(n_vox), 1, n_vox)),
    list(matrix(rnorm(n_vox), 1, n_vox), matrix(rnorm(n_vox), 1, n_vox))
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
  expect_true(file.exists(file.path(tmpdir, "beta_stability.png")))
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
    "num_passes_converged",
    "final_rho_noise_projection",
    "final_beta_stability",
    "ljung_box_p_median",
    "final_adaptive_hyperparameters"
  ) %in% names(cert)))
})
