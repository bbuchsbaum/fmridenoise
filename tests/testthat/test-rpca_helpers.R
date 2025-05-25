test_that(".ndx_rpca_single_run basic functionality", {
  td <- .generate_rpca_test_data(T_run = 20, V = 5, N_runs = 1)
  res <- NULL
  expect_no_error({
    res <- .ndx_basic_rpca_single_run(td$Y_residuals_cat, k_target = 2)
  })
  expect_true(is.list(res))
  expect_true(is.matrix(res$V_r))
  expect_equal(nrow(res$V_r), td$V)
  expect_true(ncol(res$V_r) <= 2)
  expect_true(is.matrix(res$S_matrix_TxV))
  expect_equal(dim(res$S_matrix_TxV), dim(td$Y_residuals_cat))
  expect_length(res$spike_TR_mask, td$T_run)
})

test_that(".ndx_rpca_single_run handles empty input", {
  empty <- matrix(numeric(0), 0, 0)
  res <- .ndx_basic_rpca_single_run(empty, k_target = 1)
  expect_true(is.list(res))
  expect_null(res$V_r)
  expect_equal(length(res$spike_TR_mask), 0)
})

test_that(".ndx_merge_voxel_components dimensions", {
  td <- .generate_rpca_test_data(T_run = 30, V = 6, N_runs = 2)
  V1 <- matrix(rnorm(td$V * 2), td$V, 2)
  V2 <- matrix(rnorm(td$V * 2), td$V, 2)
  Vg_concat <- .ndx_basic_merge_voxel_components(list(V1, V2), k_target = 3, strategy = "concat_svd")
  expect_true(is.matrix(Vg_concat))
  expect_equal(nrow(Vg_concat), td$V)
  expect_true(ncol(Vg_concat) <= 3)

  Vg_iter <- .ndx_basic_merge_voxel_components(list(V1, V2), k_target = 3, strategy = "iterative")
  expect_true(is.matrix(Vg_iter))
  expect_equal(nrow(Vg_iter), td$V)
  expect_true(ncol(Vg_iter) <= 3)
})
