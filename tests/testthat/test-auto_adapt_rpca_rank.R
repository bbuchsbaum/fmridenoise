context("Auto_Adapt_RPCA_Rank - k_elbow logic")

test_that("Auto_Adapt_RPCA_Rank finds elbow and respects clamps", {
  sv <- c(10, 8, 6, 0.1, 0.05)
  k <- Auto_Adapt_RPCA_Rank(sv, drop_ratio = 0.05, k_min = 2, k_max = 10)
  expect_equal(k, 4L)
})

test_that("Auto_Adapt_RPCA_Rank clamps to k_max when no drop below threshold", {
  sv <- seq(10, 1, length.out = 10)
  k <- Auto_Adapt_RPCA_Rank(sv, drop_ratio = 0.001, k_min = 2, k_max = 8)
  expect_equal(k, 8L)
})

test_that("Auto_Adapt_RPCA_Rank does not exceed length(singular_values)", {
  sv <- c(10, 5, 2)
  k <- Auto_Adapt_RPCA_Rank(sv, drop_ratio = 0.5, k_min = 5, k_max = 50)
  expect_equal(k, length(sv))
})
