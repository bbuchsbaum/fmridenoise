context("Auto_Adapt_RPCA_Rank - RPCA Rank Adaptation Logic")

test_that("Auto_Adapt_RPCA_Rank handles typical singular value spectrum", {
  singular_values1 <- 100 * exp(-(0:50)/5) # Exponential decay
  k1 <- Auto_Adapt_RPCA_Rank(singular_values1, drop_ratio = 0.02, k_min = 5, k_max = 30)
  expect_true(k1 >= 5 && k1 <= 30)
  # For this specific decay, with drop_ratio=0.02 (1/50), ratios are exp(-(0:50)/5)
  # log(0.02) = -3.91. So we need -(k-1)/5 < -3.91 => (k-1)/5 > 3.91 => k-1 > 19.55 => k > 20.55. So k=21.
  # Then clamped by k_min=5, k_max=30. Should be 21.
  expect_equal(k1, 21L)

  singular_values2 <- c(100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5, 1, 0.5)
  k2 <- Auto_Adapt_RPCA_Rank(singular_values2, drop_ratio = 0.15, k_min = 2, k_max = 10)
  # ratios: 1, .9, .8, .7, .6, .5, .4, .3, .2, .1, .05, .01, .005
  # First ratio < 0.15 is 0.1 (10th value, index 10)
  # k_candidate = 10. Clamped: max(10,2)=10, min(10,10)=10. min(10, 13)=10.
  expect_equal(k2, 10L)
})

test_that("Auto_Adapt_RPCA_Rank handles edge cases for singular_values input", {
  expect_warning(k_null <- Auto_Adapt_RPCA_Rank(NULL, k_min = 5, k_max = 10), 
                 "Auto_Adapt_RPCA_Rank: singular_values is NULL or empty.")
  expect_equal(k_null, 5L)
  
  expect_warning(k_empty <- Auto_Adapt_RPCA_Rank(numeric(0), k_min = 5, k_max = 10),
                 "Auto_Adapt_RPCA_Rank: singular_values is NULL or empty.")
  expect_equal(k_empty, 5L)

  expect_error(Auto_Adapt_RPCA_Rank(as.character(1:10)), "singular_values must be numeric")
  
  expect_warning(k_na_first <- Auto_Adapt_RPCA_Rank(c(NA, 10, 5), k_min = 3), 
                 "Auto_Adapt_RPCA_Rank: first singular value is not finite.")
  expect_equal(k_na_first, 3L)
  
  expect_warning(k_inf_first <- Auto_Adapt_RPCA_Rank(c(Inf, 10, 5), k_min = 3),
                 "Auto_Adapt_RPCA_Rank: first singular value is not finite.")
  expect_equal(k_inf_first, 3L)
})

test_that("Auto_Adapt_RPCA_Rank respects k_min and k_max clamping", {
  singular_values <- c(100, 10, 1, 0.1) # Ratios: 1, 0.1, 0.01, 0.001
  
  # drop_ratio = 0.05 -> k_candidate should be 3 (where ratio is 0.01)
  expect_equal(Auto_Adapt_RPCA_Rank(singular_values, drop_ratio = 0.05, k_min = 1, k_max = 4), 3L)
  
  # k_min forces higher
  expect_equal(Auto_Adapt_RPCA_Rank(singular_values, drop_ratio = 0.05, k_min = 4, k_max = 4), 4L) 
  expect_equal(Auto_Adapt_RPCA_Rank(singular_values, drop_ratio = 0.05, k_min = 5, k_max = 10), 4L) # Clamped by length(sv)

  # k_max forces lower
  expect_equal(Auto_Adapt_RPCA_Rank(singular_values, drop_ratio = 0.0001, k_min = 1, k_max = 2), 2L) # drop_ratio would pick 4, but k_max=2
  expect_equal(Auto_Adapt_RPCA_Rank(singular_values, drop_ratio = 0.5, k_min = 1, k_max = 1), 1L)   # drop_ratio would pick 2, but k_max=1
})

test_that("Auto_Adapt_RPCA_Rank handles cases where drop_ratio is never/always met", {
  singular_values <- c(100, 90, 80, 70) # Ratios all > 0.5
  # Never meets drop_ratio = 0.1, should return length(singular_values) or k_max
  expect_equal(Auto_Adapt_RPCA_Rank(singular_values, drop_ratio = 0.1, k_min = 1, k_max = 10), length(singular_values)) # 4
  expect_equal(Auto_Adapt_RPCA_Rank(singular_values, drop_ratio = 0.1, k_min = 1, k_max = 3), 3L) # clamped by k_max

  # Always meets drop_ratio = 0.95 (except first one)
  # k_candidate will be 2. Then clamped by k_min.
  expect_equal(Auto_Adapt_RPCA_Rank(singular_values, drop_ratio = 0.95, k_min = 1, k_max = 10), 2L)
  expect_equal(Auto_Adapt_RPCA_Rank(singular_values, drop_ratio = 0.95, k_min = 3, k_max = 10), 3L)
})

test_that("Auto_Adapt_RPCA_Rank result is always integer and within overall bounds", {
  sv <- 100 * exp(-(0:60)/8)
  k_test <- Auto_Adapt_RPCA_Rank(sv, drop_ratio = 0.01, k_min = 10, k_max = 40)
  expect_true(is.integer(k_test))
  expect_true(k_test >= 10L)
  expect_true(k_test <= 40L)
  expect_true(k_test <= length(sv))
}) 