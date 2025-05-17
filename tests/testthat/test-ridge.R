context("ndx_solve_ridge - Functionality")

# Helper function to generate some data
.generate_ridge_test_data <- function(n_timepoints, n_regressors, n_voxels, true_betas, noise_sd = 1.0, na_head = 0) {
  X <- matrix(rnorm(n_timepoints * n_regressors), n_timepoints, n_regressors)
  # Ensure true_betas is a matrix of correct dimensions (n_regressors x n_voxels)
  if (is.vector(true_betas)) {
    if (length(true_betas) != n_regressors) stop("Length of true_betas vector must match n_regressors")
    true_betas_matrix <- matrix(rep(true_betas, n_voxels), nrow = n_regressors, ncol = n_voxels)
  } else if (is.matrix(true_betas)) {
    if (nrow(true_betas) != n_regressors || ncol(true_betas) != n_voxels) {
      stop("Dimensions of true_betas matrix are incorrect.")
    }
    true_betas_matrix <- true_betas
  } else {
    stop("true_betas must be a vector or a matrix.")
  }
  
  Y_signal <- X %*% true_betas_matrix
  Y_noise <- matrix(rnorm(n_timepoints * n_voxels, sd = noise_sd), n_timepoints, n_voxels)
  Y <- Y_signal + Y_noise
  
  na_mask_vec <- rep(FALSE, n_timepoints)
  if (na_head > 0 && n_timepoints > na_head) {
    Y[1:na_head, ] <- NA
    # X[1:na_head, ] <- NA # Not strictly needed to make X have NAs for this test, Y is enough
    na_mask_vec[1:na_head] <- TRUE
  }
  
  return(list(Y = Y, X = X, true_betas = true_betas_matrix, na_mask = na_mask_vec))
}

test_that("ndx_solve_ridge computes correct betas for a simple case", {
  set.seed(123)
  n_tp <- 100
  n_reg <- 2
  n_vox <- 1
  true_b <- c(2, -3)
  lambda <- 1.0
  
  test_data <- .generate_ridge_test_data(n_tp, n_reg, n_vox, true_b, noise_sd = 0.1)
  
  # Manual calculation for single voxel
  # beta = (X'X + lambda * I)^(-1) X'Y
  XtX <- crossprod(test_data$X)
  I <- diag(n_reg)
  expected_betas_manual <- solve(XtX + lambda * I) %*% crossprod(test_data$X, test_data$Y[,1,drop=FALSE])
  
  ridge_betas <- ndx_solve_ridge(test_data$Y, test_data$X, lambda_ridge = lambda)
  
  expect_true(is.matrix(ridge_betas))
  expect_equal(nrow(ridge_betas), n_reg)
  expect_equal(ncol(ridge_betas), n_vox)
  expect_equal(ridge_betas[,1], as.vector(expected_betas_manual), tolerance = 1e-6)
})

test_that("ndx_solve_ridge handles lambda_ridge = 0 (OLS-like)", {
  set.seed(456)
  n_tp <- 50
  n_reg <- 3
  n_vox <- 2
  true_b_mat <- matrix(c(1,0.5,-1, -2,1,0.3), nrow=n_reg, ncol=n_vox) 
  
  test_data <- .generate_ridge_test_data(n_tp, n_reg, n_vox, true_b_mat, noise_sd = 0.2)
  
  # OLS solution: (X'X)^-1 X'Y
  # Using lm.fit for robustness against collinearity, though unlikely here
  ols_betas_list <- lapply(1:n_vox, function(v_idx) coef(lm.fit(test_data$X, test_data$Y[,v_idx])))
  ols_betas_manual <- do.call(cbind, ols_betas_list)
  dimnames(ols_betas_manual) <- NULL # Strip dimnames for comparison consistency

  ridge_betas_lambda0 <- ndx_solve_ridge(test_data$Y, test_data$X, lambda_ridge = 0)
  
  expect_equal(ridge_betas_lambda0, ols_betas_manual, tolerance = 1e-6)
})

test_that("ndx_solve_ridge shrinks betas with high lambda", {
  set.seed(789)
  n_tp <- 50
  n_reg <- 2
  n_vox <- 1
  true_b <- c(10, -15)
  test_data <- .generate_ridge_test_data(n_tp, n_reg, n_vox, true_b, noise_sd = 1)

  ridge_betas_low_lambda <- ndx_solve_ridge(test_data$Y, test_data$X, lambda_ridge = 0.1)
  ridge_betas_high_lambda <- ndx_solve_ridge(test_data$Y, test_data$X, lambda_ridge = 10000)
  
  # Check that sum of absolute beta values is much smaller for high lambda
  expect_lt(sum(abs(ridge_betas_high_lambda)), sum(abs(ridge_betas_low_lambda)) * 0.1)
  # Check that high lambda betas are close to zero
  expect_true(all(abs(ridge_betas_high_lambda) < 0.1))
})

test_that("ndx_solve_ridge handles na_mask and NAs in Y_whitened correctly", {
  set.seed(101)
  n_tp <- 60
  n_reg <- 2
  n_vox <- 2
  true_b_mat <- matrix(c(1,2, -1,0.5), n_reg, n_vox)
  na_count = 3
  
  test_data_na <- .generate_ridge_test_data(n_tp, n_reg, n_vox, true_b_mat, noise_sd = 0.1, na_head = na_count)
  Y_with_na <- test_data_na$Y
  X_for_na <- test_data_na$X
  na_mask_provided <- test_data_na$na_mask
  
  # 1. With na_mask provided
  betas_with_mask <- ndx_solve_ridge(Y_with_na, X_for_na, lambda_ridge = 1.0, na_mask = na_mask_provided)
  
  # Manual calculation on clean data
  Y_clean_manual <- Y_with_na[!na_mask_provided, , drop = FALSE]
  X_clean_manual <- X_for_na[!na_mask_provided, , drop = FALSE]
  XtX_clean <- crossprod(X_clean_manual)
  I_clean <- diag(n_reg)
  expected_betas_clean_manual <- solve(XtX_clean + 1.0 * I_clean) %*% crossprod(X_clean_manual, Y_clean_manual)
  
  expect_equal(betas_with_mask, expected_betas_clean_manual, tolerance = 1e-6)

  # 2. Without na_mask (should auto-detect and remove NA rows from Y)
  betas_auto_na <- ndx_solve_ridge(Y_with_na, X_for_na, lambda_ridge = 1.0, na_mask = NULL)
  expect_equal(betas_auto_na, expected_betas_clean_manual, tolerance = 1e-6, 
               label = "Should match manual clean when na_mask is NULL and Y has NAs")
               
  # 3. No NAs present, na_mask=NULL (should run on full data)
  test_data_no_na <- .generate_ridge_test_data(n_tp, n_reg, n_vox, true_b_mat, noise_sd = 0.1, na_head = 0)
  XtX_full <- crossprod(test_data_no_na$X)
  expected_betas_full_manual <- solve(XtX_full + 1.0 * diag(n_reg)) %*% crossprod(test_data_no_na$X, test_data_no_na$Y)
  betas_no_na_null_mask <- ndx_solve_ridge(test_data_no_na$Y, test_data_no_na$X, lambda_ridge = 1.0, na_mask = NULL)
  expect_equal(betas_no_na_null_mask, expected_betas_full_manual, tolerance = 1e-6)
  
  # 4. All Y rows are NA
  Y_all_na <- matrix(NA_real_, nrow=n_tp, ncol=n_vox)
  expect_warning(
    betas_all_na <- ndx_solve_ridge(Y_all_na, X_for_na, lambda_ridge = 1.0, na_mask = NULL),
    "No timepoints remaining after NA removal. Returning NA betas."
  )
  expect_true(all(is.na(betas_all_na)))
  expect_equal(nrow(betas_all_na), n_reg)
  expect_equal(ncol(betas_all_na), n_vox)
})

test_that("ndx_solve_ridge input validation works", {
  Y_good <- matrix(rnorm(20*2), 20, 2)
  X_good <- matrix(rnorm(20*3), 20, 3)
  lambda_good <- 1.0
  na_mask_good <- rep(FALSE, 20)
  
  expect_error(ndx_solve_ridge(as.data.frame(Y_good), X_good, lambda_good), "Y_whitened must be a numeric matrix.")
  expect_error(ndx_solve_ridge(Y_good, as.data.frame(X_good), lambda_good), "X_whitened must be a numeric matrix.")
  expect_error(ndx_solve_ridge(Y_good, X_good[1:10,], lambda_good), "Y_whitened and X_whitened must have the same number of rows")
  expect_error(ndx_solve_ridge(Y_good, X_good, "a"), "lambda_ridge must be a single non-negative numeric value.")
  expect_error(ndx_solve_ridge(Y_good, X_good, -1), "lambda_ridge must be a single non-negative numeric value.")
  expect_error(ndx_solve_ridge(Y_good, X_good, c(1,2)), "lambda_ridge must be a single non-negative numeric value.")
  expect_error(ndx_solve_ridge(Y_good, X_good, lambda_good, na_mask = rep(FALSE, 10)), "Provided na_mask must be a logical vector of length equal to nrows of Y_whitened.")
  expect_error(ndx_solve_ridge(Y_good, X_good, lambda_good, na_mask = as.numeric(na_mask_good)), "Provided na_mask must be a logical vector of length equal to nrows of Y_whitened.")

  # More regressors than timepoints after NA removal
  X_too_many_reg <- matrix(rnorm(5*6), 5, 6)
  Y_short <- matrix(rnorm(5*1), 5, 1)
  expect_warning(
    b_warn <- ndx_solve_ridge(Y_short, X_too_many_reg, lambda_good, na_mask=rep(FALSE,5)),
    "Number of timepoints after NA removal \\(5\\) is less than number of regressors \\(6\\)"
  )
  expect_true(all(is.na(b_warn)))
  
  # Zero regressors or voxels
  expect_warning(b_zero_reg <- ndx_solve_ridge(Y_good, X_good[,0,drop=FALSE], lambda_good), "X_whitened has 0 regressors.")
  expect_equal(dim(b_zero_reg), c(0, ncol(Y_good)))
  expect_warning(b_zero_vox <- ndx_solve_ridge(Y_good[,0,drop=FALSE], X_good, lambda_good), "Y_whitened has 0 voxels.")
  expect_equal(dim(b_zero_vox), c(ncol(X_good), 0))
})

context("ndx_extract_task_betas - Functionality")

test_that("ndx_extract_task_betas extracts specified betas correctly", {
  betas_full <- matrix(1:12, nrow = 4, ncol = 3)
  colnames_x <- paste0("reg", 1:4) # Matches rows of betas_full if rownames were set
  rownames(betas_full) <- colnames_x # Simulate betas from a named X
  
  tasks_to_extract <- c("reg2", "reg4")
  extracted <- ndx_extract_task_betas(betas_full, colnames_x, tasks_to_extract)
  
  expect_true(!is.null(extracted))
  expect_equal(nrow(extracted), length(tasks_to_extract))
  expect_equal(ncol(extracted), ncol(betas_full))
  expect_equal(rownames(extracted), tasks_to_extract)
  expect_equal(extracted["reg2", ], betas_full["reg2", ])
  expect_equal(extracted["reg4", ], betas_full["reg4", ])
  
  # Test with X_whitened_colnames not matching rownames of betas_whitened (should still work using X_whitened_colnames for indexing)
  betas_mismatched_rownames <- matrix(1:12, nrow = 4, ncol = 3)
  rownames(betas_mismatched_rownames) <- paste0("original_beta_name", 1:4) # Give it different rownames
  colnames_x <- paste0("reg", 1:4) # Target names for X
  tasks_to_extract <- c("reg2", "reg4")

  expect_warning(
    extracted_mismatch_rn <- ndx_extract_task_betas(betas_mismatched_rownames, colnames_x, tasks_to_extract),
    "Rownames of betas_whitened do not match X_whitened_colnames"
  )
  expect_true(!is.null(extracted_mismatch_rn))
  # We expect indexing to be based on colnames_x, so reg2 is the 2nd element, reg4 is the 4th
  expect_equal(extracted_mismatch_rn["reg2", ], betas_mismatched_rownames[2, ]) 
  expect_equal(extracted_mismatch_rn["reg4", ], betas_mismatched_rownames[4, ]) 
})

test_that("ndx_extract_task_betas handles missing task regressors", {
  betas_full <- matrix(1:8, nrow = 4, ncol = 2)
  colnames_x <- c("taskA", "noise1", "taskB", "noise2")
  rownames(betas_full) <- colnames_x
  
  # One missing, one present
  tasks_mixed <- c("taskA", "taskC_missing")
  expect_warning(
    extracted_mixed <- ndx_extract_task_betas(betas_full, colnames_x, tasks_mixed),
    "The following task regressors were not found in X_whitened_colnames: taskC_missing"
  )
  expect_true(!is.null(extracted_mixed))
  expect_equal(nrow(extracted_mixed), 1)
  expect_equal(rownames(extracted_mixed), "taskA")
  expect_equal(extracted_mixed["taskA", ], betas_full["taskA", ])
  
  # All missing
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
  
  # Zero row betas
  betas_zero_row <- matrix(numeric(0), ncol=2, nrow=0) # Ensure it's a numeric matrix
  colnames_zero_row <- character(0)
  expect_warning(ndx_extract_task_betas(betas_zero_row, colnames_zero_row, tasks_good), "Input betas_whitened has 0 regressors.")
  expect_null(suppressWarnings(ndx_extract_task_betas(betas_zero_row, colnames_zero_row, tasks_good)))
}) 