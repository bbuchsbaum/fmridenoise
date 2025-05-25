context("ndx_progressive_viz functions")

test_that("ndx_generate_progressive_enhancement creates widget and returns iframe", {
  tmpdir <- tempfile("progressive")
  dir.create(tmpdir)

  workflow_mock <- list(
    annihilation_mode_active = TRUE,
    gdlite_pcs = matrix(rnorm(10), nrow = 5),
    rpca_orthogonalized = matrix(rnorm(10), nrow = 5),
    spectral_orthogonalized = matrix(rnorm(10), nrow = 5),
    diagnostics_per_pass = list(list(DES = 0.7))
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

  workflow_mock <- list(
    annihilation_mode_active = TRUE,
    gdlite_pcs = matrix(rnorm(10), nrow = 5),
    rpca_orthogonalized = matrix(rnorm(10), nrow = 5),
    spectral_orthogonalized = matrix(rnorm(10), nrow = 5),
    diagnostics_per_pass = list(list(DES = 0.7))
  )

  html_lines <- c("<html>", "<body>", "<p>Report</p>", "</body>", "</html>")

  updated_html <- ndx_add_progressive_enhancement_to_report(html_lines, workflow_mock, output_dir = tmpdir)

  expect_true(file.exists(file.path(tmpdir, "progressive_enhancement.html")))
  expect_true(any(grepl("Progressive Enhancement Visualization", updated_html)))
  expect_gt(length(updated_html), length(html_lines))
})
