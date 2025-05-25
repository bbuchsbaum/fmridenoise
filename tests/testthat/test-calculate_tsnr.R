context("calculate_tsnr")

test_that("basic tSNR computation", {
  set.seed(123)
  noise <- matrix(rnorm(20 * 2, sd = 2), nrow = 20, ncol = 2)
  signal <- matrix(5, nrow = 20, ncol = 2)
  Y <- signal + noise
  tsnr <- calculate_tsnr(Y, detrend = FALSE)
  manual <- colMeans(Y) / apply(Y, 2, sd)
  expect_equal(tsnr, manual, tolerance = 1e-6)
})

test_that("detrending removes linear trends", {
  time <- 1:10
  Y <- cbind(time, rnorm(10))
  tsnr_dt <- calculate_tsnr(Y, detrend = TRUE)
  expect_true(is.na(tsnr_dt[1]))
  expect_true(!is.na(tsnr_dt[2]))
})

test_that("zero variance voxels produce NA", {
  Y <- matrix(3, nrow = 5, ncol = 2)
  tsnr <- calculate_tsnr(Y)
  expect_true(all(is.na(tsnr)))
})
