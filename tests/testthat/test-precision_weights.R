context("ndx_precision_weights_from_S")

test_that("precision weights computed correctly", {
  S <- matrix(c(-2, 0, 2, 4), nrow = 2)
  w <- ndx_precision_weights_from_S(S, mad_floor = 1e-6)
  mad_S <- max(mad(abs(S), na.rm = TRUE), 1e-6)
  expected <- 1 / (abs(S) / mad_S + 1)^2
  expect_equal(w, expected)
  expect_equal(dim(w), dim(S))
})

test_that("NA and Inf values yield weight of one", {
  S <- matrix(c(1, NA, Inf, -2), nrow = 2)
  w <- ndx_precision_weights_from_S(S, mad_floor = 1e-6, na_zero = TRUE)
  abs_S <- abs(S)
  abs_S[!is.finite(abs_S)] <- 0
  mad_S <- max(mad(abs_S, na.rm = TRUE), 1e-6)
  expected <- 1 / (abs_S / mad_S + 1)^2
  expected[!is.finite(expected)] <- 1
  expect_equal(w, expected)
  expect_true(all(is.finite(w)))
})
