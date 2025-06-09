context("AR filter helper functions")

# Helper to manually filter a vector with given AR coefficients
.manual_ar_filter <- function(x, coeffs) {
  p <- length(coeffs)
  n <- length(x)
  res <- rep(NA_real_, n)
  for (t in seq_len(n)) {
    if (t > p) {
      val <- x[t]
      for (k in seq_len(p)) {
        val <- val - coeffs[k] * x[t - k]
      }
      res[t] <- val
    }
  }
  res
}


test_that(".apply_ar_filter_to_matrix_cols filters correctly and handles short matrices", {
  M <- matrix(1:12, nrow = 6, ncol = 2)
  coeff <- c(0.5, -0.25)
  res <- ndx:::.apply_ar_filter_to_matrix_cols(M, coeff, ar_order = 2L)
  expected <- cbind(
    .manual_ar_filter(M[,1], coeff),
    .manual_ar_filter(M[,2], coeff)
  )
  expect_equal(res, expected)

  short_M <- matrix(rnorm(4), nrow = 2, ncol = 2)
  expect_warning(res_short <- ndx:::.apply_ar_filter_to_matrix_cols(short_M, coeff, ar_order = 2L),
                 "Matrix to be filtered")
  expect_true(all(is.na(res_short)))
})


test_that(".apply_ar_filter_voxelwise applies coefficients per voxel", {
  set.seed(1)
  Y <- matrix(rnorm(30), nrow = 10, ncol = 3)
  coeffs <- rbind(c(0.6, -0.2), c(0, 0), c(0.3, 0.1))
  res <- ndx:::.apply_ar_filter_voxelwise(Y, coeffs, ar_order = 2L)
  expected <- cbind(
    .manual_ar_filter(Y[,1], coeffs[1,]),
    Y[,2],
    .manual_ar_filter(Y[,3], coeffs[3,])
  )
  expect_equal(res, expected)
})


test_that(".fit_ar_single_voxel estimates coefficients and handles edge cases", {
  set.seed(42)
  phi_true <- c(0.6, -0.2)
  series <- as.numeric(arima.sim(n = 200, list(ar = phi_true)))
  fit <- ndx:::.fit_ar_single_voxel(series, order = 2L)
  fit_ref <- stats::ar.yw(series, aic = FALSE, order.max = 2)
  expect_equal(fit[1:2], fit_ref$ar, tolerance = 1e-6)
  expect_equal(fit[3], fit_ref$var.pred, tolerance = 1e-6)

  const_series <- rep(1, 50)
  fit_const <- ndx:::.fit_ar_single_voxel(const_series, order = 2L)
  expect_equal(fit_const[1:2], c(0,0))
  expect_true(is.na(fit_const[3]))

  weights_zero <- rep(0, length(series))
  fit_zero_w <- ndx:::.fit_ar_single_voxel(series, order = 2L, weights_vec = weights_zero)
  expect_equal(fit_zero_w[1:2], c(0,0))
  expect_true(is.na(fit_zero_w[3]))
})

