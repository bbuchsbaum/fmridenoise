context("ndx_progressive_viz functions")

test_that("ndx_generate_progressive_enhancement creates widget and returns iframe", {
  tmpdir <- tempfile("progressive")
  dir.create(tmpdir)

  # Create mock data with all required fields
  n_timepoints <- 50
  n_voxels <- 10
  n_regressors <- 5
  
  workflow_mock <- list(
    annihilation_mode_active = TRUE,
    gdlite_pcs = matrix(rnorm(n_timepoints * 3), nrow = n_timepoints, ncol = 3),
    rpca_orthogonalized = matrix(rnorm(n_timepoints * 2), nrow = n_timepoints, ncol = 2),
    spectral_orthogonalized = matrix(rnorm(n_timepoints * 2), nrow = n_timepoints, ncol = 2),
    diagnostics_per_pass = list(
      list(DES = 0.5),
      list(DES = 0.7),
      list(DES = 0.8)
    ),
    beta_history_per_pass = list(
      matrix(rnorm(n_regressors * n_voxels), nrow = n_regressors, ncol = n_voxels),
      matrix(rnorm(n_regressors * n_voxels), nrow = n_regressors, ncol = n_voxels),
      matrix(rnorm(n_regressors * n_voxels), nrow = n_regressors, ncol = n_voxels)
    ),
    pass0_residuals = matrix(rnorm(n_timepoints * n_voxels), nrow = n_timepoints, ncol = n_voxels),
    Y_residuals_final_unwhitened = matrix(rnorm(n_timepoints * n_voxels), nrow = n_timepoints, ncol = n_voxels)
  )

  iframe_html <- NULL
  expect_no_error({
    iframe_html <- ndx_generate_progressive_enhancement(workflow_mock, output_dir = tmpdir)
  })

  expect_true(file.exists(file.path(tmpdir, "progressive_enhancement.html")))
  expect_true(grepl("<iframe", iframe_html, fixed = TRUE))
})

test_that("ndx_add_progressive_enhancement_to_report inserts section", {
  tmpdir <- tempfile("progressive")
  dir.create(tmpdir)

  # Create mock data with all required fields
  n_timepoints <- 50
  n_voxels <- 10
  n_regressors <- 5
  
  workflow_mock <- list(
    annihilation_mode_active = TRUE,
    gdlite_pcs = matrix(rnorm(n_timepoints * 3), nrow = n_timepoints, ncol = 3),
    rpca_orthogonalized = matrix(rnorm(n_timepoints * 2), nrow = n_timepoints, ncol = 2),
    spectral_orthogonalized = matrix(rnorm(n_timepoints * 2), nrow = n_timepoints, ncol = 2),
    diagnostics_per_pass = list(
      list(DES = 0.5),
      list(DES = 0.7),
      list(DES = 0.8)
    ),
    beta_history_per_pass = list(
      matrix(rnorm(n_regressors * n_voxels), nrow = n_regressors, ncol = n_voxels),
      matrix(rnorm(n_regressors * n_voxels), nrow = n_regressors, ncol = n_voxels),
      matrix(rnorm(n_regressors * n_voxels), nrow = n_regressors, ncol = n_voxels)
    ),
    pass0_residuals = matrix(rnorm(n_timepoints * n_voxels), nrow = n_timepoints, ncol = n_voxels),
    Y_residuals_final_unwhitened = matrix(rnorm(n_timepoints * n_voxels), nrow = n_timepoints, ncol = n_voxels)
  )

  html_lines <- c("<html>", "<body>", "<p>Report</p>", "</body>", "</html>")

  updated_html <- ndx_add_progressive_enhancement_to_report(html_lines, workflow_mock, output_dir = tmpdir)

  expect_true(file.exists(file.path(tmpdir, "progressive_enhancement.html")))
  expect_true(any(grepl("Progressive Enhancement Visualization", updated_html)))
  expect_gt(length(updated_html), length(html_lines))
})
