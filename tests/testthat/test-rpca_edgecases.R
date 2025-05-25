context("ndx_rpca_temporal_components_multirun - Edge Cases")

# This file adds additional tests to ensure the RPCA workflow
# is robust to unusual or extreme inputs.

# k_global_target = 0 should simply return NULL

test_that("k_global_target of zero returns NULL", {
  td <- .generate_rpca_test_data(N_runs = 1)
  res <- ndx_rpca_temporal_components_multirun(
    Y_residuals_cat = td$Y_residuals_cat,
    run_idx = td$run_idx,
    k_global_target = 0,
    user_options = list()
  )
  expect_null(res)
})

# Negative k_global_target should error

test_that("negative k_global_target errors", {
  td <- .generate_rpca_test_data(N_runs = 1)
  expect_error(
    ndx_rpca_temporal_components_multirun(
      Y_residuals_cat = td$Y_residuals_cat,
      run_idx = td$run_idx,
      k_global_target = -1,
      user_options = list()
    ),
    "k_global_target must be non-negative"
  )
})

# k_per_run_target <= 0 should default to k_global_target

test_that("non-positive k_per_run_target defaults to k_global_target", {
  td <- .generate_rpca_test_data(N_runs = 1)
  opts <- list(k_per_run_target = 0)
  res <- ndx_rpca_temporal_components_multirun(
    Y_residuals_cat = td$Y_residuals_cat,
    run_idx = td$run_idx,
    k_global_target = 2,
    user_options = opts
  )
  expect_true(!is.null(res))
  expect_equal(ncol(res$C_components), 2)
  expect_length(res$spike_TR_mask, td$total_T)
})

# Constant residual matrix should not cause errors

test_that("constant residual matrix is handled", {
  td <- .generate_rpca_test_data(N_runs = 1)
  const_mat <- matrix(3, nrow = td$total_T, ncol = td$V)
  res <- NULL
  expect_no_error({
    res <- ndx_rpca_temporal_components_multirun(
      Y_residuals_cat = const_mat,
      run_idx = td$run_idx,
      k_global_target = 1,
      user_options = list()
    )
  })
  expect_true(!is.null(res))
  expect_true(is.matrix(res$C_components))
  expect_equal(nrow(res$C_components), td$total_T)
  expect_true(ncol(res$C_components) <= 1)
})
