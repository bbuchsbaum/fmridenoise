context("ndx_extract_gdlite_pcs")

test_that("valid extraction with and without mask", {
  set.seed(123)
  Y <- matrix(rnorm(20 * 10), nrow = 20, ncol = 10)
  X <- matrix(1, nrow = 20, ncol = 1)

  pcs_nomask <- ndx_extract_gdlite_pcs(Y, X, n_pcs = 2)
  expect_true(is.matrix(pcs_nomask))
  expect_equal(dim(pcs_nomask), c(20, 2))

  mask <- rep(c(TRUE, FALSE), length.out = 10)
  pcs_mask <- ndx_extract_gdlite_pcs(Y, X, n_pcs = 2, voxel_mask = mask)
  expect_true(is.matrix(pcs_mask))
  expect_equal(dim(pcs_mask), c(20, 2))
})

test_that("empty noise pool due to mask returns NULL", {
  set.seed(1)
  Y <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5)
  X <- matrix(1, nrow = 10, ncol = 1)
  mask_none <- rep(FALSE, 5)
  expect_warning(res <- ndx_extract_gdlite_pcs(Y, X, n_pcs = 2, voxel_mask = mask_none))
  expect_null(res)
})

test_that("zero variance voxels return NULL", {
  Y <- matrix(5, nrow = 8, ncol = 4)  # constant data
  X <- matrix(1, nrow = 8, ncol = 1)
  expect_warning(res <- ndx_extract_gdlite_pcs(Y, X, n_pcs = 1))
  expect_null(res)
})

test_that("invalid inputs trigger errors", {
  Y <- matrix(rnorm(10), nrow = 5)
  X <- matrix(1, nrow = 5, ncol = 1)
  expect_error(ndx_extract_gdlite_pcs(list(1,2,3), X, n_pcs = 2))
  expect_error(ndx_extract_gdlite_pcs(Y, list(), n_pcs = 2))
})
