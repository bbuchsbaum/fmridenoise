context("workflow helper utilities")

test_that("ndx_validate_process_subject_inputs works", {
  Y <- matrix(0, nrow = 5, ncol = 2)
  events <- data.frame(onsets = c(1), durations = c(1), condition = factor("A"), blockids = 1)
  motion <- matrix(0, nrow = 5, ncol = 1)
  expect_silent(ndx_validate_process_subject_inputs(Y, events, motion, rep(1,5), 1))
  expect_error(ndx_validate_process_subject_inputs(1, events, motion, rep(1,5), 1))
  expect_error(ndx_validate_process_subject_inputs(Y, events[0,], motion, rep(1,5), 1))
})

test_that("ndx_prepare_workflow_options merges defaults", {
  user_opts <- list(max_passes = 5, opts_pass0 = list(poly_degree = 3))
  opts <- ndx_prepare_workflow_options(user_opts)
  expect_equal(opts$max_passes, 5)
  expect_equal(opts$opts_pass0$poly_degree, 3)
  expect_true(!is.null(opts$opts_hrf))
})

test_that("ndx_run_annihilation_setup handles disabled mode", {
  Y <- matrix(rnorm(20), nrow = 10, ncol = 2)
  events <- data.frame(onsets = 1, durations = 1, condition = factor("A"), blockids = 1)
  motion <- matrix(0, nrow = 10, ncol = 1)
  res <- ndx_run_annihilation_setup(Y, events, motion, rep(1,10), 1, list(annihilation_enable_mode = FALSE))
  expect_null(res$selected_pcs)
})

test_that("ndx_run_annihilation_setup returns PCs and diagnostics when enabled", {
  set.seed(123)
  n_time <- 6
  run_idx <- rep(1:2, each = 3)
  Y <- matrix(rnorm(n_time * 4), nrow = n_time, ncol = 4)
  events <- data.frame(
    onsets = c(0, 3),
    durations = c(1, 1),
    condition = factor("A"),
    blockids = 1:2
  )
  motion <- matrix(rnorm(n_time * 2), nrow = n_time, ncol = 2)
  opts <- list(
    annihilation_enable_mode = TRUE,
    annihilation_gdlite_poly_degree = 0,
    annihilation_gdlite_k_max = 2,
    annihilation_gdlite_r2_thresh_noise_pool = 1,
    annihilation_gdlite_tsnr_thresh_noise_pool = -Inf,
    annihilation_gdlite_r2_thresh_good_voxels = -Inf
  )

  res <- ndx_run_annihilation_setup(Y, events, motion, run_idx, 1, opts, verbose = FALSE)

  expect_true(is.matrix(res$selected_pcs))
  expect_true(ncol(res$selected_pcs) > 0)
  expect_true("r2_vals_by_k" %in% names(res$gdlite_results))

  opts$annihilation_enable_mode <- FALSE
  res_dis <- ndx_run_annihilation_setup(Y, events, motion, run_idx, 1, opts, verbose = FALSE)
  expect_null(res_dis$selected_pcs)
})
