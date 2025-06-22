context("additional utility functions")

# Test calculate_beta_stability

test_that("calculate_beta_stability computes pairwise correlations", {
  betas_pass1 <- list(matrix(1:4, nrow = 2), matrix(1:4, nrow = 2))
  betas_pass2 <- list(matrix(c(1,0,0,1), nrow = 2), matrix(c(-1,0,0,-1), nrow = 2))
  res <- calculate_beta_stability(list(betas_pass1, betas_pass2))
  expect_equal(res[1], 1, tolerance = 1e-8)
  expect_equal(res[2], -1, tolerance = 1e-8)
})

test_that("calculate_beta_stability returns NA with single run", {
  betas_single <- list(list(matrix(1, nrow = 1, ncol = 1)))
  res <- calculate_beta_stability(betas_single)
  expect_true(is.na(res[1]))
})

# Tests for Ljung-Box helpers

test_that("compute_ljung_box_pvalues returns p-values per series", {
  set.seed(1)
  Y <- matrix(rnorm(100 * 2), nrow = 100, ncol = 2)
  pvals <- compute_ljung_box_pvalues(Y, lag = 5)
  expect_length(pvals, 2)
  expect_true(all(pvals > 0 & pvals <= 1))
})

test_that("ndx_ljung_box_pval handles matrix and vector inputs", {
  set.seed(2)
  Y <- matrix(rnorm(60 * 3), nrow = 60, ncol = 3)
  p_matrix <- ndx_ljung_box_pval(Y, lag = 5)
  p_vector <- ndx_ljung_box_pval(Y[,1], lag = 5)
  expect_true(is.numeric(p_matrix) && length(p_matrix) == 1)
  expect_true(is.numeric(p_vector) && length(p_vector) == 1)
  expect_true(p_matrix > 0 && p_matrix <= 1)
  expect_true(p_vector > 0 && p_vector <= 1)
})

test_that("ljung box helpers return NA for constant data", {
  Y_const <- matrix(1, nrow = 10, ncol = 2)
  pvals <- compute_ljung_box_pvalues(Y_const, lag = 5)
  expect_true(all(is.na(pvals)))
  expect_true(is.na(ndx_ljung_box_pval(Y_const, lag = 5)))
})

test_that("ljung box helpers handle NAs and short series", {
  set.seed(3)
  Y <- matrix(rnorm(8 * 2), nrow = 8, ncol = 2)
  Y[1:4, 2] <- NA  # too few non-NA observations for column 2
  pvals <- compute_ljung_box_pvalues(Y, lag = 5)
  expect_true(!is.na(pvals[1]) && is.na(pvals[2]))

  short_vec <- rnorm(4)
  expect_true(is.na(ndx_ljung_box_pval(short_vec, lag = 5)))
})

test_that("merge_lists recursively merges defaults and user opts", {
  defaults <- list(a = 1, b = 2, c = list(d = 3, e = 4))
  user <- list(b = 20, c = list(e = 40))
  res <- merge_lists(defaults, user)
  expect_equal(res$a, 1)
  expect_equal(res$b, 20)
  expect_equal(res$c$d, 3)
  expect_equal(res$c$e, 40)

  res2 <- merge_lists(defaults, NULL)
  expect_equal(res2, defaults)
})

test_that("calculate_beta_stability validates input list", {
  expect_error(calculate_beta_stability(list()), "betas_per_pass must be a non-empty list")
  expect_error(calculate_beta_stability(matrix(1:4,2,2)), "betas_per_pass must be a non-empty list")
})

# Test ndx_default_user_options structure

test_that("ndx_default_user_options returns expected structure", {
  opts <- ndx_default_user_options()
  expect_true(is.list(opts))
  expect_true(all(c("max_passes", "opts_pass0", "opts_hrf", "opts_ridge", "opts_annihilation") %in% names(opts)))
  expect_true(is.numeric(opts$max_passes) && opts$max_passes > 0)
  expect_true(is.list(opts$opts_hrf) && is.list(opts$opts_ridge))
})

