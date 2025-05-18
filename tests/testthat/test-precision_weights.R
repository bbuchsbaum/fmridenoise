context("ndx_precision_weights_from_S")

test_that("precision weights computed correctly", {
  S <- matrix(c(-2, 0, 2, 4), nrow = 2)
  w <- ndx_precision_weights_from_S(S, mad_floor = 1e-6)
  mad_S <- max(mad(abs(S), na.rm = TRUE), 1e-6)
  expected <- 1 / (abs(S) / mad_S + 1)^2
  expect_equal(w, expected)
  expect_equal(dim(w), dim(S))
})
