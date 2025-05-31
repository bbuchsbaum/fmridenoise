context("ridge input helper functions")

# Use helper to generate simple data
set.seed(42)
data <- .create_ridge_test_data(10, 2, 1, c(1, -1), noise_sd = 0.1)
Y <- data$Y
X <- data$X


# ---- validate_ridge_inputs ----

test_that("validate_ridge_inputs detects invalid inputs and accepts valid ones", {
  K_diag <- rep(0.1, ncol(X))
  expect_true(ndx:::validate_ridge_inputs(Y, X, K_penalty_diag = K_diag))

  # mismatched rows
  expect_warning(
    res_mismatch <- ndx:::validate_ridge_inputs(Y[-1,,drop=FALSE], X, K_diag),
    "same number of rows"
  )
  expect_false(res_mismatch)

  # negative penalty values
  expect_warning(
    res_neg <- ndx:::validate_ridge_inputs(Y, X, K_penalty_diag = c(-0.1, 0.1)),
    "non-negative"
  )
  expect_false(res_neg)

  # weights vector wrong length
  expect_warning(
    res_wlen <- ndx:::validate_ridge_inputs(Y, X, K_penalty_diag = K_diag,
                                            weights = rep(1, nrow(Y)-1)),
    "weights vector length"
  )
  expect_false(res_wlen)

  # use_penalty_matrix without matrices
  expect_warning(
    res_pm <- ndx:::validate_ridge_inputs(Y, X, K_penalty_diag = NULL,
                                          projection_mats = NULL,
                                          lambda_values = NULL,
                                          use_penalty_matrix = TRUE),
    "Either K_penalty_mat or projection_mats"
  )
  expect_false(res_pm)
})

# ---- apply_na_mask_and_weights ----

test_that("apply_na_mask_and_weights removes rows and subsets weights", {
  na_mask <- rep(FALSE, nrow(Y))
  na_mask[1:2] <- TRUE
  w <- matrix(seq_len(nrow(Y)), ncol = 1)

  res <- ndx:::apply_na_mask_and_weights(Y, X, na_mask = na_mask, weights = w)
  expect_equal(nrow(res$Y), nrow(Y) - 2)
  expect_equal(nrow(res$X), nrow(X) - 2)
  expect_equal(nrow(res$W), nrow(Y) - 2)
  expect_equal(res$W[,1], w[!na_mask,1])

  # all masked returns NULL
  na_all <- rep(TRUE, nrow(Y))
  expect_warning(res_all <- ndx:::apply_na_mask_and_weights(Y, X, na_mask = na_all, weights = w),
                 "All timepoints are masked")
  expect_null(res_all)

  # invalid mask length -> uses all rows
  expect_warning(res_bad <- ndx:::apply_na_mask_and_weights(Y, X, na_mask = na_mask[-1], weights = w),
                 "na_mask must be")
  expect_equal(nrow(res_bad$Y), nrow(Y))
  expect_equal(nrow(res_bad$W), nrow(Y))
})

