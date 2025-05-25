context("ndx_solve_anisotropic_ridge - Functionality")
# Using .create_ridge_test_data from helper-ridge_data.R


test_that("ndx_solve_anisotropic_ridge computes correct betas for a simple case", {
  set.seed(123)
  n_tp <- 100
  n_reg <- 2
  n_vox <- 1
  true_b <- c(2, -3)
  lambda_val <- 0.5
  
  test_data <- .create_ridge_test_data(n_tp, n_reg, n_vox, true_b, noise_sd = 0.1)
  K_diag <- rep(lambda_val, ncol(test_data$X))
  
  betas_ridge <- ndx_solve_anisotropic_ridge(test_data$Y, test_data$X, K_penalty_diag = K_diag)$betas
  dimnames(betas_ridge) <- NULL # Strip dimnames for comparison

  XtX <- crossprod(test_data$X)
  K_matrix_manual <- diag(lambda_val, ncol(test_data$X))
  expected_betas_manual <- solve(XtX + K_matrix_manual) %*% crossprod(test_data$X, test_data$Y[,1,drop=FALSE])
  dimnames(expected_betas_manual) <- NULL # Strip dimnames
  
  expect_true(!is.null(betas_ridge))
  expect_equal(dim(betas_ridge), c(n_reg, n_vox))
  expect_equal(as.vector(betas_ridge[,1]), as.vector(expected_betas_manual), tolerance = 1e-6)
})

test_that("ndx_solve_anisotropic_ridge works with full penalty matrix", {
  set.seed(123)
  n_tp <- 100
  n_reg <- 2
  n_vox <- 1
  true_b <- c(2, -3)
  lambda_val <- 0.5

  td <- .create_ridge_test_data(n_tp, n_reg, n_vox, true_b, noise_sd = 0.1)
  K_mat <- diag(lambda_val, n_reg)

  betas <- ndx_solve_anisotropic_ridge(td$Y, td$X,
                                       K_penalty_mat = K_mat,
                                       use_penalty_matrix = TRUE)$betas

  expected <- solve(crossprod(td$X) + K_mat) %*% crossprod(td$X, td$Y)
  dimnames(betas) <- NULL
  dimnames(expected) <- NULL
  expect_equal(betas, expected, tolerance = 1e-6)
})

test_that("ndx_solve_anisotropic_ridge handles K_penalty_diag = 0 (OLS-like)", {
  set.seed(456)
  n_tp <- 50
  n_reg <- 3
  n_vox <- 2
  true_b_mat <- matrix(c(1,0.5,-1, -2,1,0.3), nrow=n_reg, ncol=n_vox) 
  
  test_data <- .create_ridge_test_data(n_tp, n_reg, n_vox, true_b_mat, noise_sd = 0.2)
  K_diag_zero <- rep(0, ncol(test_data$X))
  
  betas_ols_like <- ndx_solve_anisotropic_ridge(test_data$Y, test_data$X, K_penalty_diag = K_diag_zero)$betas
  dimnames(betas_ols_like) <- NULL # Strip dimnames
  
  ols_betas_lm_fit <- stats::lm.fit(test_data$X, test_data$Y)$coefficients
  dimnames(ols_betas_lm_fit) <- NULL 
  
  expect_equal(betas_ols_like, ols_betas_lm_fit, tolerance = 1e-6)
})

test_that("ndx_solve_anisotropic_ridge shrinks betas with high K_penalty_diag values", {
  set.seed(789)
  n_tp <- 50
  n_reg <- 2
  n_vox <- 1
  true_b <- c(10, -15)
  test_data <- .create_ridge_test_data(n_tp, n_reg, n_vox, true_b, noise_sd = 1)

  K_diag_low <- rep(0.1, ncol(test_data$X))
  K_diag_high <- rep(10000, ncol(test_data$X))

  betas_low_lambda <- ndx_solve_anisotropic_ridge(test_data$Y, test_data$X, K_penalty_diag = K_diag_low)$betas
  betas_high_lambda <- ndx_solve_anisotropic_ridge(test_data$Y, test_data$X, K_penalty_diag = K_diag_high)$betas
  
  expect_lt(sum(abs(betas_high_lambda)), sum(abs(betas_low_lambda)) * 0.1)
  expect_true(all(abs(betas_high_lambda) < 0.1))
})

test_that("ndx_solve_anisotropic_ridge handles na_mask correctly", {
  set.seed(101)
  n_tp <- 60
  n_reg <- 2
  n_vox <- 2
  true_b_mat <- matrix(c(1,2, -1,0.5), n_reg, n_vox)
  na_count = 3
  
  test_data_na <- .create_ridge_test_data(n_tp, n_reg, n_vox, true_b_mat, noise_sd = 0.1, na_head = na_count)
  Y_with_na <- test_data_na$Y
  X_for_na <- test_data_na$X
  na_mask_provided <- test_data_na$na_mask
  K_diag <- rep(1.0, ncol(X_for_na))
  
  betas_with_mask <- ndx_solve_anisotropic_ridge(Y_with_na, X_for_na, K_penalty_diag = K_diag, na_mask = na_mask_provided)$betas
  
  Y_clean_manual <- Y_with_na[!na_mask_provided, , drop = FALSE]
  X_clean_manual <- X_for_na[!na_mask_provided, , drop = FALSE]
  K_matrix_manual <- diag(1.0, ncol(X_clean_manual))
  expected_betas_clean_manual <- solve(crossprod(X_clean_manual) + K_matrix_manual, 
                                       crossprod(X_clean_manual, Y_clean_manual))
  
  expect_equal(betas_with_mask, expected_betas_clean_manual, tolerance = 1e-6)

  betas_raw_na <- ndx_solve_anisotropic_ridge(Y_with_na, X_for_na, K_penalty_diag = K_diag, na_mask = NULL)$betas
  expect_true(any(is.na(betas_raw_na))) 

  all_na_mask <- rep(TRUE, nrow(Y_with_na))
  expect_warning(
    betas3 <- ndx_solve_anisotropic_ridge(Y_with_na, X_for_na, K_penalty_diag = K_diag, na_mask = all_na_mask),
    "All timepoints are masked by na_mask"
  )
  expect_null(betas3)
})


test_that("ndx_solve_anisotropic_ridge input validation works", {
  Y_good <- matrix(rnorm(20*2), 20, 2)
  X_good <- matrix(rnorm(20*3), 20, 3)
  K_diag_good <- rep(0.1, ncol(X_good))
  
  expect_warning(ndx_solve_anisotropic_ridge(as.data.frame(Y_good), X_good, K_diag_good), "Y_whitened must be a numeric matrix", fixed = TRUE)
  expect_warning(ndx_solve_anisotropic_ridge(Y_good, as.data.frame(X_good), K_diag_good), "X_whitened must be a numeric matrix", fixed = TRUE)
  expect_warning(ndx_solve_anisotropic_ridge(Y_good, X_good[1:10, ], K_diag_good), "Y_whitened and X_whitened must have the same number of rows", fixed = TRUE)
  expect_warning(ndx_solve_anisotropic_ridge(Y_good, X_good, K_diag_good[1]), "K_penalty_diag must be a numeric vector with length equal to ncol(X_whitened).", fixed = TRUE)
  expect_warning(ndx_solve_anisotropic_ridge(Y_good, X_good, "a"), "K_penalty_diag must be a numeric vector", fixed = TRUE) 
  expect_warning(ndx_solve_anisotropic_ridge(Y_good, X_good, c(-0.1, 0.1, 0.1)), "All elements of K_penalty_diag must be non-negative", fixed = TRUE)
  expect_warning(ndx_solve_anisotropic_ridge(Y_good, X_good, K_diag_good, na_mask=rep(TRUE, 5)), "na_mask must be a logical vector with length equal to nrow(Y_whitened).", fixed = TRUE)

  w_bad_neg <- c(rep(1, nrow(Y_good) - 1), -1)
  expect_warning(res_neg <- ndx_solve_anisotropic_ridge(Y_good, X_good, K_diag_good, weights = w_bad_neg),
                 "weights must contain non-negative finite numbers", fixed = TRUE)
  expect_null(res_neg)

  w_bad_inf <- c(rep(1, nrow(Y_good) - 1), Inf)
  expect_warning(res_inf <- ndx_solve_anisotropic_ridge(Y_good, X_good, K_diag_good, weights = w_bad_inf),
                 "weights must contain non-negative finite numbers", fixed = TRUE)
  expect_null(res_inf)
  
  Y_short_for_warn <- Y_good[1:2, , drop=FALSE]
  X_wide_for_warn <- X_good[1:2, , drop=FALSE] 
  expect_warning(ndx_solve_anisotropic_ridge(Y_short_for_warn, X_wide_for_warn, K_diag_good), 
                 "Number of timepoints after NA removal (2) is less than number of regressors (3)", fixed = TRUE)
  
  expect_warning(ndx_solve_anisotropic_ridge(Y_good, X_good[,0,drop=FALSE], K_diag_good), "K_penalty_diag must be a numeric vector with length equal to ncol(X_whitened).", fixed = TRUE) 
  expect_null(ndx_solve_anisotropic_ridge(Y_good, X_good[,0,drop=FALSE], numeric(0)))
})

context("ndx_extract_task_betas - Functionality")

test_that("ndx_extract_task_betas extracts specified betas correctly", {
  betas_full <- matrix(1:12, nrow = 4, ncol = 3)
  colnames_x <- paste0("reg", 1:4) 
  rownames(betas_full) <- colnames_x 
  
  tasks_to_extract <- c("reg2", "reg4")
  extracted <- ndx_extract_task_betas(betas_full, colnames_x, tasks_to_extract)
  
  expect_true(!is.null(extracted))
  expect_equal(nrow(extracted), length(tasks_to_extract))
  expect_equal(ncol(extracted), ncol(betas_full))
  expect_equal(rownames(extracted), tasks_to_extract)
  expect_equal(extracted["reg2", ], betas_full["reg2", ])
  expect_equal(extracted["reg4", ], betas_full["reg4", ])
  
  betas_mismatched_rownames <- matrix(1:12, nrow = 4, ncol = 3)
  rownames(betas_mismatched_rownames) <- paste0("original_beta_name", 1:4)
  
  expect_warning(
    extracted_mismatch_rn <- ndx_extract_task_betas(betas_mismatched_rownames, colnames_x, tasks_to_extract),
    "Rownames of betas_whitened do not match X_whitened_colnames"
  )
  expect_true(!is.null(extracted_mismatch_rn))
  expect_equal(extracted_mismatch_rn["reg2", ], betas_mismatched_rownames[2, ]) 
  expect_equal(extracted_mismatch_rn["reg4", ], betas_mismatched_rownames[4, ]) 
})

test_that("ndx_extract_task_betas handles missing task regressors", {
  betas_full <- matrix(1:8, nrow = 4, ncol = 2)
  colnames_x <- c("taskA", "noise1", "taskB", "noise2")
  rownames(betas_full) <- colnames_x
  
  tasks_mixed <- c("taskA", "taskC_missing")
  expect_warning(
    extracted_mixed <- ndx_extract_task_betas(betas_full, colnames_x, tasks_mixed),
    "The following task regressors were not found in X_whitened_colnames: taskC_missing"
  )
  expect_true(!is.null(extracted_mixed))
  expect_equal(nrow(extracted_mixed), 1)
  expect_equal(rownames(extracted_mixed), "taskA")
  expect_equal(extracted_mixed["taskA", ], betas_full["taskA", ])
  
  tasks_all_missing <- c("taskX", "taskY")
  expect_warning(
    extracted_all_missing <- ndx_extract_task_betas(betas_full, colnames_x, tasks_all_missing),
    "No specified task regressors were found in X_whitened_colnames. Returning NULL."
  )
  expect_null(extracted_all_missing)
})

test_that("ndx_extract_task_betas input validation works", {
  betas_good <- matrix(1:6, 3, 2)
  colnames_good <- c("r1", "r2", "r3")
  tasks_good <- c("r1")
  
  expect_error(ndx_extract_task_betas(as.data.frame(betas_good), colnames_good, tasks_good), "betas_whitened must be a numeric matrix.")
  expect_error(ndx_extract_task_betas(betas_good, colnames_good[-1], tasks_good), "Number of rows in betas_whitened must match the length of X_whitened_colnames.")
  expect_error(ndx_extract_task_betas(betas_good, as.list(colnames_good), tasks_good), "X_whitened_colnames must be a character vector.")
  expect_error(ndx_extract_task_betas(betas_good, colnames_good, character(0)), "task_regressor_names must be a non-empty character vector.")
  expect_error(ndx_extract_task_betas(betas_good, colnames_good, list("r1")), "task_regressor_names must be a non-empty character vector.")
  
  betas_zero_row <- matrix(numeric(0), ncol=2, nrow=0) 
  colnames_zero_row <- character(0)
  expect_warning(ndx_extract_task_betas(betas_zero_row, colnames_zero_row, tasks_good), "Input betas_whitened has 0 regressors.")
  expect_null(suppressWarnings(ndx_extract_task_betas(betas_zero_row, colnames_zero_row, tasks_good)))
})

test_that("ndx_compute_projection_matrices constructs valid projectors", {
  set.seed(1)
  U_noise <- matrix(rnorm(20), nrow = 4)
  proj <- ndx_compute_projection_matrices(U_Noise = U_noise, n_regressors = 4)
  expect_true(all(dim(proj$P_Noise) == c(4,4)))
  expect_equal(proj$P_Noise + proj$P_Signal, diag(4), tolerance = 1e-6)
  expect_true(max(abs(proj$P_Noise %*% proj$P_Signal)) < 1e-6)
})

test_that("ndx_update_lambda_aggressiveness adjusts lambda", {
  lam <- 1
  expect_gt(ndx_update_lambda_aggressiveness(lam, 0.2), lam)
  expect_lt(ndx_update_lambda_aggressiveness(lam, 0.0), lam)
})

test_that("ndx_solve_anisotropic_ridge works with projection matrices", {
  td <- .create_ridge_test_data(50, 2, 1, c(2, -1), noise_sd = 0.1)
  U_noise <- matrix(rnorm(100), nrow = 2)
  proj <- ndx_compute_projection_matrices(U_Noise = U_noise, n_regressors = 2)
  lambda_vals <- list(lambda_parallel = 0.5, lambda_perp_signal = 0.05)
  betas <- ndx_solve_anisotropic_ridge(td$Y, td$X, K_penalty_diag = NULL,
                                       projection_mats = proj,
                                       lambda_values = lambda_vals)$betas
  expect_equal(dim(betas), c(2,1))
})

test_that("full penalty from projection matrices matches manual", {
  td <- .create_ridge_test_data(30, 2, 1, c(2, -1), noise_sd = 0.1)
  U_noise <- diag(2)[,1,drop=FALSE]
  proj <- ndx_compute_projection_matrices(U_Noise = U_noise, n_regressors = 2)
  lambda_vals <- list(lambda_parallel = 0.5, lambda_perp_signal = 0.05)

  K_manual <- lambda_vals$lambda_parallel * proj$P_Noise +
              lambda_vals$lambda_perp_signal * proj$P_Signal

  betas_pm <- ndx_solve_anisotropic_ridge(td$Y, td$X,
                                          projection_mats = proj,
                                          lambda_values = lambda_vals,
                                          use_penalty_matrix = TRUE)$betas

  expected <- solve(crossprod(td$X) + K_manual) %*% crossprod(td$X, td$Y)
  dimnames(betas_pm) <- NULL
  dimnames(expected) <- NULL
  expect_equal(betas_pm, expected, tolerance = 1e-6)
})

test_that("ndx_solve_anisotropic_ridge handles weights", {
  data <- .create_ridge_test_data(40, 2, 1, c(1, -1), noise_sd = 0.1)
  K_diag <- rep(0.1, 2)
  w <- matrix(1, nrow = nrow(data$Y), ncol = 1)
  betas_unw <- ndx_solve_anisotropic_ridge(data$Y, data$X, K_penalty_diag = K_diag)$betas
  betas_w <- ndx_solve_anisotropic_ridge(data$Y, data$X, K_penalty_diag = K_diag, weights = w)$betas
  expect_equal(betas_unw, betas_w)
})
