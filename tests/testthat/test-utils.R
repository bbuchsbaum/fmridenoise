context("Utility functions from R/ndx_utils.R")

test_that("calculate_residuals_ols computes correct residuals", {
  Y <- matrix(c(1,2,3, 2,4,6, 3,6,9), ncol=3, byrow=TRUE)
  X <- matrix(c(1,1, 1,2, 1,3), ncol=2, byrow=TRUE)
  
  # Manually calculate residuals for a simple case
  # Y = X %*% beta + E
  # beta_hat = solve(t(X) %*% X) %*% t(X) %*% Y
  # residuals = Y - X %*% beta_hat
  
  # For first column of Y: y1 = c(1,2,3)
  # beta1_hat = solve(t(X)%*%X) %*% t(X) %*% y1
  XtX <- t(X) %*% X
  XtY1 <- t(X) %*% Y[,1,drop=FALSE]
  beta1_hat <- solve(XtX, XtY1)
  residuals1_manual <- Y[,1,drop=FALSE] - X %*% beta1_hat
  
  # For second column of Y: y2 = c(2,4,6) (which is 2*y1)
  # beta2_hat should be 2*beta1_hat
  # residuals2_manual should be 2*residuals1_manual
  XtY2 <- t(X) %*% Y[,2,drop=FALSE]
  beta2_hat <- solve(XtX, XtY2)
  residuals2_manual <- Y[,2,drop=FALSE] - X %*% beta2_hat
  
  # For third column of Y: y3 = c(3,6,9) (which is 3*y1)
  XtY3 <- t(X) %*% Y[,3,drop=FALSE]
  beta3_hat <- solve(XtX, XtY3)
  residuals3_manual <- Y[,3,drop=FALSE] - X %*% beta3_hat
  
  expected_residuals <- cbind(residuals1_manual, residuals2_manual, residuals3_manual)
  colnames(expected_residuals) <- NULL # lm.fit residuals don't have names by default
  rownames(expected_residuals) <- NULL # lm.fit residuals don't have names by default

  # Test with calculate_residuals_ols
  res_ols <- calculate_residuals_ols(Y, X)
  
  expect_equal(res_ols, expected_residuals, tolerance = 1e-8)
  expect_true(is.matrix(res_ols))
  expect_equal(dim(res_ols), dim(Y))
})

test_that("calculate_residuals_ols handles rank-deficient X and gives warning", {
  Y <- matrix(rnorm(10), ncol=2)
  X_rankdef <- matrix(c(1,1, 1,1, 1,1, 1,1, 1,1), ncol=2) # Rank 1, 2 cols
  
  expect_warning(
    res_ols_rankdef <- calculate_residuals_ols(Y, X_rankdef),
    "Design matrix is rank deficient. Rank = 1 Columns = 2"
  )
  
  # Residuals should still be calculable
  # Use lm.fit directly to get expected residuals, mirroring calculate_residuals_ols internal
  expected_fit_rankdef <- stats::lm.fit(X_rankdef, Y)
  expected_residuals_rankdef <- expected_fit_rankdef$residuals
  # Input Y has no dimnames, so lm.fit$residuals should not have them either.
  # If there's still an attribute issue, res_ols_rankdef might need dimnames(.) <- NULL
  
  expect_equal(res_ols_rankdef, expected_residuals_rankdef, tolerance = 1e-8)
})

test_that("calculate_residuals_ols handles X with zero columns and gives warning", {
  Y <- matrix(rnorm(10), ncol=2)
  X_zero_cols <- matrix(numeric(0), nrow=5, ncol=0)
  
  expect_warning(
    res_ols_zerocols <- calculate_residuals_ols(Y, X_zero_cols),
    "Design matrix X has zero columns. Returning Y as residuals."
  )
  expect_equal(res_ols_zerocols, Y)
})

test_that("calculate_residuals_ols errors with mismatched rows in Y and X", {
  Y <- matrix(rnorm(10), ncol=2) # 5 rows
  X <- matrix(rnorm(12), ncol=2) # 6 rows
  expect_error(
    calculate_residuals_ols(Y, X),
    "Number of rows in Y and X must match for OLS."
  )
})

test_that("calculate_R2_voxelwise computes correct R2 values", {
  Y_obs <- matrix(c(1,2,3,4, 2,4,6,8, 10,20,30,40), ncol=3, byrow=FALSE)
  # Col1: linear trend, perfect fit expected if model captures it.
  # Col2: 2*Col1
  # Col3: 10*Col1
  
  # Case 1: Perfect fit (residuals are zero)
  Y_res_perfect <- matrix(0, nrow=4, ncol=3)
  expected_R2_perfect <- c(1, 1, 1)
  expect_equal(calculate_R2_voxelwise(Y_obs, Y_res_perfect), expected_R2_perfect, tolerance = 1e-9)
  
  # Case 2: No fit (residuals are Y_obs - mean(Y_obs) if model is just intercept, or Y_obs if model is null)
  # If Y_residuals are Y_obs (model explains nothing beyond what Y_obs itself is, relative to 0)
  # and Y_obs has non-zero variance, R2 should be < 1. If Y_obs is centered, R2 would be 0.
  # calculate_R2_voxelwise assumes Y_residuals are from a model fit to Y_observed.
  # TSS is based on Y_observed. If RSS = TSS (i.e., model explains nothing beyond mean), R2=0.
  # If model makes things worse (RSS > TSS), R2 < 0, which is then clamped to 0.

  # Test with residuals = Y_obs (implies model predicted zero everywhere)
  # This should lead to R2 = 1 - sum(Y_obs^2) / sum((Y_obs - mean(Y_obs))^2)
  # This isn't standard R2, usually residuals are Y_obs - Y_hat.
  # Let's create more meaningful residuals.
  
  # Model Y_hat = mean(Y_obs) for each voxel
  Y_hat_mean <- apply(Y_obs, 2, function(col) rep(mean(col), length(col)))
  Y_res_mean <- Y_obs - Y_hat_mean
  expected_R2_mean <- c(0,0,0) # Model only explains the mean, so R2=0 conventionally.
  expect_equal(calculate_R2_voxelwise(Y_obs, Y_res_mean), expected_R2_mean, tolerance = 1e-9)
  
  # Case 3: Partial fit
  # For Y_obs[,1] = c(1,2,3,4), mean=2.5, TSS = (1-2.5)^2+(2-2.5)^2+(3-2.5)^2+(4-2.5)^2 = 2.25+0.25+0.25+2.25 = 5
  # Say Y_hat is c(1.5, 2, 2.5, 3)
  Y_hat_partial_col1 <- c(1.5, 2, 2.5, 3)
  Y_res_partial_col1 <- Y_obs[,1] - Y_hat_partial_col1 # c(-0.5, 0, 0.5, 1)
  RSS_col1 <- sum(Y_res_partial_col1^2) # 0.25 + 0 + 0.25 + 1 = 1.5
  R2_col1 <- 1 - (RSS_col1 / 5) # 1 - 1.5/5 = 1 - 0.3 = 0.7
  
  Y_res_partial <- cbind(Y_res_partial_col1, 2*Y_res_partial_col1, 10*Y_res_partial_col1)
  expected_R2_partial <- c(R2_col1, R2_col1, R2_col1) # R2 is scale invariant
  expect_equal(calculate_R2_voxelwise(Y_obs, Y_res_partial), expected_R2_partial, tolerance = 1e-9)
  
  # Case 4: Residuals make R2 negative (clamped to 0)
  Y_res_worse <- Y_obs * 2 # Residuals are much larger than original data variance from mean
  expected_R2_worse <- c(0,0,0)
  expect_equal(calculate_R2_voxelwise(Y_obs, Y_res_worse), expected_R2_worse, tolerance = 1e-9)
  
  # Case 5: Y_observed has zero variance (TSS = 0)
  Y_obs_flat <- matrix(rep(c(2,4,20), each=4), ncol=3, byrow=FALSE) # Each col is constant
  Y_res_flat <- matrix(0, nrow=4, ncol=3)
  expected_R2_flat <- c(0,0,0) # TSS is 0, R2 should be 0
  expect_equal(calculate_R2_voxelwise(Y_obs_flat, Y_res_flat), expected_R2_flat, tolerance = 1e-9)

  Y_res_flat_noise <- matrix(rnorm(12), ncol=3) # Non-zero residuals
  expect_equal(calculate_R2_voxelwise(Y_obs_flat, Y_res_flat_noise), expected_R2_flat, tolerance = 1e-9)
})

test_that("calculate_R2_voxelwise handles NA values correctly", {
  Y_obs_na <- matrix(c(1,2,NA,4, 2,NA,6,8, 10,20,30,NA), ncol=3, byrow=FALSE)
  Y_res_na <- matrix(rnorm(12), ncol=3) # Assume residuals are complete for simplicity of R2 calc
  
  # For col1: Y_obs=c(1,2,NA,4), Y_res=Y_res_na[,1]. n_obs=3. mean_obs=(1+2+4)/3 = 7/3.
  # TSS1 = (1-7/3)^2 + (2-7/3)^2 + (4-7/3)^2 = (-4/3)^2 + (-1/3)^2 + (5/3)^2 = (16+1+25)/9 = 42/9 = 14/3
  # RSS1 = sum(Y_res_na[c(1,2,4),1]^2)
  # R2_1 = 1 - RSS1/TSS1
  
  # Manual calculation for first column with NAs
  y_obs1_valid <- Y_obs_na[!is.na(Y_obs_na[,1]), 1]
  y_res1_valid <- Y_res_na[!is.na(Y_obs_na[,1]), 1] # Align residuals with valid observed data
  tss1_manual <- sum((y_obs1_valid - mean(y_obs1_valid))^2)
  rss1_manual <- sum(y_res1_valid^2)
  r2_1_manual <- if (abs(tss1_manual) < 1e-9) 0 else max(0, 1 - rss1_manual / tss1_manual)
  
  # Manual calculation for second column with NAs
  y_obs2_valid <- Y_obs_na[!is.na(Y_obs_na[,2]), 2]
  y_res2_valid <- Y_res_na[!is.na(Y_obs_na[,2]), 2]
  tss2_manual <- sum((y_obs2_valid - mean(y_obs2_valid))^2)
  rss2_manual <- sum(y_res2_valid^2)
  r2_2_manual <- if (abs(tss2_manual) < 1e-9) 0 else max(0, 1 - rss2_manual / tss2_manual)
  
  # Manual calculation for third column with NAs
  y_obs3_valid <- Y_obs_na[!is.na(Y_obs_na[,3]), 3]
  y_res3_valid <- Y_res_na[!is.na(Y_obs_na[,3]), 3]
  tss3_manual <- sum((y_obs3_valid - mean(y_obs3_valid))^2)
  rss3_manual <- sum(y_res3_valid^2)
  r2_3_manual <- if (abs(tss3_manual) < 1e-9) 0 else max(0, 1 - rss3_manual / tss3_manual)
  
  expected_R2_na <- c(r2_1_manual, r2_2_manual, r2_3_manual)
  expect_equal(calculate_R2_voxelwise(Y_obs_na, Y_res_na), expected_R2_na, tolerance = 1e-9)
  
  # Test where residuals also have NAs (should align with Y_obs_na for RSS calculation)
  Y_res_na_aligned <- Y_res_na
  Y_res_na_aligned[is.na(Y_obs_na)] <- NA # Make residuals NA where Y_obs is NA
  
  rss1_manual_aligned <- sum(Y_res_na_aligned[!is.na(Y_obs_na[,1]), 1]^2, na.rm=TRUE)
  r2_1_manual_aligned <- if (abs(tss1_manual) < 1e-9) 0 else max(0, 1 - rss1_manual_aligned / tss1_manual)
  rss2_manual_aligned <- sum(Y_res_na_aligned[!is.na(Y_obs_na[,2]), 2]^2, na.rm=TRUE)
  r2_2_manual_aligned <- if (abs(tss2_manual) < 1e-9) 0 else max(0, 1 - rss2_manual_aligned / tss2_manual)
  rss3_manual_aligned <- sum(Y_res_na_aligned[!is.na(Y_obs_na[,3]), 3]^2, na.rm=TRUE)
  r2_3_manual_aligned <- if (abs(tss3_manual) < 1e-9) 0 else max(0, 1 - rss3_manual_aligned / tss3_manual)
  
  expected_R2_na_aligned <- c(r2_1_manual_aligned, r2_2_manual_aligned, r2_3_manual_aligned)
  expect_equal(calculate_R2_voxelwise(Y_obs_na, Y_res_na_aligned), expected_R2_na_aligned, tolerance = 1e-9)
})

test_that("calculate_R2_voxelwise errors on mismatched dimensions", {
  Y_obs <- matrix(rnorm(10), ncol=2)
  Y_res_wrong_rows <- matrix(rnorm(8), ncol=2)
  Y_res_wrong_cols <- matrix(rnorm(10), ncol=1)
  expect_error(calculate_R2_voxelwise(Y_obs, Y_res_wrong_rows), "Dimensions of Y_observed and Y_residuals must match.")
  expect_error(calculate_R2_voxelwise(Y_obs, Y_res_wrong_cols), "Dimensions of Y_observed and Y_residuals must match.")
}) 
test_that("calculate_R2_voxelwise handles all-NA columns and extra NA residuals", {
  Y_obs <- matrix(c(1, NA, 3, 4,  NA, NA, NA, NA), ncol = 2)
  Y_res <- matrix(c(0.1, 0.2, NA, -0.1,  0.5, -0.5, 0.2, -0.2), ncol = 2)

  # Manual R2 for first column ignoring NA residual
  mask1 <- !is.na(Y_obs[,1])
  y_obs1 <- Y_obs[mask1, 1]
  y_res1 <- Y_res[mask1, 1]
  tss1 <- sum((y_obs1 - mean(y_obs1))^2)
  rss1 <- sum(y_res1^2, na.rm = TRUE)
  expected_r2_col1 <- if (abs(tss1) < 1e-9) 0 else max(0, 1 - rss1 / tss1)

  r2_vals <- calculate_R2_voxelwise(Y_obs, Y_res)
  expect_equal(r2_vals[1], expected_r2_col1, tolerance = 1e-9)
  expect_equal(r2_vals[2], 0)
})
