context("Global AR stability safeguard")

test_that("Unstable averaged AR coefficients skip design whitening", {
  set.seed(42)
  n_tp <- 60
  Y <- matrix(rnorm(n_tp * 2), nrow = n_tp, ncol = 2)
  X <- matrix(rnorm(n_tp * 2), nrow = n_tp, ncol = 2)
  Y_res <- Y
  fake_coeffs <- c(1.2, 0.8)  # purposely unstable

  stub_colMeans <- function(x, na.rm = FALSE, dims = 1) {
    if (is.matrix(x) && identical(dim(x), c(2L, 2L))) {
      return(fake_coeffs)
    } else {
      base::colMeans(x, na.rm = na.rm, dims = dims)
    }
  }

  ndx_env <- new.env(parent = environment(ndx_ar2_whitening))
  ndx_env$colMeans <- stub_colMeans
  ndx_ar2_whitening_stub <- ndx_ar2_whitening
  environment(ndx_ar2_whitening_stub) <- ndx_env

  expect_warning(
    res <- ndx_ar2_whitening_stub(Y, X, Y_res, order = 2L, verbose = FALSE),
    regexp = "Averaged AR coefficients"
  )
  expect_null(res$AR_coeffs_global)
  expect_identical(res$X_whitened, X)
})
