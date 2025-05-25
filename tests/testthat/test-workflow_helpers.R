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
