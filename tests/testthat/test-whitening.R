context("ndx_ar2_whitening - Functionality")

# Helper function to generate AR(order) data
.generate_ar_data <- function(n_timepoints, ar_coeffs, innov_sd = 1.0, burn_in = 100) {
  order <- length(ar_coeffs)
  total_points <- n_timepoints + burn_in
  innovations <- rnorm(total_points, mean = 0, sd = innov_sd)
  series <- numeric(total_points)
  
  for (i in (order + 1):total_points) {
    ar_part <- sum(ar_coeffs * series[(i - 1):(i - order)])
    series[i] <- ar_part + innovations[i]
  }
  return(list(series = series[(burn_in + 1):total_points], innovations = innovations[(burn_in + 1):total_points]))
}

test_that("ndx_ar2_whitening estimates AR(2) coeffs and whitens Y_data correctly", {
  set.seed(123)
  n_voxels <- 2
  n_timepoints <- 200
  ar_order <- 2

  # Voxel 1: Known AR(2) process
  true_ar1_coeffs <- c(0.6, -0.3)
  innov_sd1 <- 1.5
  ar_data1 <- .generate_ar_data(n_timepoints, true_ar1_coeffs, innov_sd1)
  
  # Voxel 2: Different AR(2) process
  true_ar2_coeffs <- c(0.4, 0.2)
  innov_sd2 <- 0.8
  ar_data2 <- .generate_ar_data(n_timepoints, true_ar2_coeffs, innov_sd2)

  Y_data_sim <- cbind(ar_data1$series, ar_data2$series)
  # For AR fitting, use the data itself if it's already supposed to be residuals
  # Or use actual residuals from a model fit if Y_data contains signal + noise
  Y_residuals_for_AR_fit_sim <- Y_data_sim 

  # Dummy design matrix (will be whitened by global AR)
  X_design_sim <- matrix(rnorm(n_timepoints * 3), n_timepoints, 3)
  colnames(X_design_sim) <- paste0("reg", 1:3)

  whitening_results <- NULL
  expect_no_error({
    whitening_results <- ndx_ar2_whitening(Y_data_sim, X_design_sim, Y_residuals_for_AR_fit_sim, order = ar_order)
  })
  
  expect_true(is.list(whitening_results))
  expect_named(whitening_results, c("Y_whitened", "X_whitened", "AR_coeffs_voxelwise", "AR_coeffs_global", "var_innovations_voxelwise", "na_mask"))

  # Check AR_coeffs_voxelwise
  coeffs_est <- whitening_results$AR_coeffs_voxelwise
  expect_equal(nrow(coeffs_est), n_voxels)
  expect_equal(ncol(coeffs_est), ar_order)
  # Check if estimated coefficients are reasonably close to true ones
  expect_true(all(abs(coeffs_est[1,] - true_ar1_coeffs) < 0.15), label = "Voxel 1 AR coeffs close to true")
  expect_true(all(abs(coeffs_est[2,] - true_ar2_coeffs) < 0.15), label = "Voxel 2 AR coeffs close to true")

  # Check var_innovations_voxelwise
  var_innov_est <- whitening_results$var_innovations_voxelwise
  expect_length(var_innov_est, n_voxels)
  expect_true(abs(var_innov_est[1] - innov_sd1^2) < innov_sd1^2 * 0.5, label = "Voxel 1 innov var close to true") # Wider tolerance for variance
  expect_true(abs(var_innov_est[2] - innov_sd2^2) < innov_sd2^2 * 0.5, label = "Voxel 2 innov var close to true")
  
  # Check Y_whitened
  Y_w <- whitening_results$Y_whitened
  expect_equal(dim(Y_w), dim(Y_data_sim))
  expect_true(all(is.na(Y_w[1:ar_order, ])), label = "First 'order' rows of Y_whitened should be NA")
  
  # Whitened series should approximate the original innovations (scaled)
  # Check acf of whitened series (should be low for lags > 0)
  na_mask_from_output <- whitening_results$na_mask # or simply use 1:ar_order
  
  # Voxel 1 ACF check
  if(any(whitening_results$AR_coeffs_voxelwise[1, ] != 0)) {
    acf_y_w1 <- acf(Y_w[!na_mask_from_output, 1], plot = FALSE, lag.max = 5)
    expect_true(all(abs(acf_y_w1$acf[-1]) < 0.25), label = "ACF of whitened Y_data (Voxel 1) should be low if coeffs are non-zero")
  } else {
    # message for testthat output if skipped
    skip("ACF check for Voxel 1 skipped as AR coefficients were zeroed out.")
  }
  
  # Voxel 2 ACF check
  if(any(whitening_results$AR_coeffs_voxelwise[2, ] != 0)) {
    acf_y_w2 <- acf(Y_w[!na_mask_from_output, 2], plot = FALSE, lag.max = 5)
    expect_true(all(abs(acf_y_w2$acf[-1]) < 0.25), label = "ACF of whitened Y_data (Voxel 2) should be low if coeffs are non-zero")
  } else {
    skip("ACF check for Voxel 2 skipped as AR coefficients were zeroed out.")
  }
  
  # Check X_whitened (whitened by global AR model)
  X_w <- whitening_results$X_whitened
  expect_equal(dim(X_w), dim(X_design_sim))
  if (!is.null(whitening_results$AR_coeffs_global) && any(whitening_results$AR_coeffs_global != 0)) {
    expect_true(!is.null(whitening_results$AR_coeffs_global), label = "Global AR coeffs for X should exist if used.")
    expect_length(whitening_results$AR_coeffs_global, ar_order)
    expect_true(all(is.na(X_w[1:ar_order, ])), label = "First 'order' rows of X_whitened should be NA if whitened")
    # Check if X_w is different from X_design_sim (unless global AR coeffs are ~0)
    expect_false(isTRUE(all.equal(X_w[!na_mask_from_output, ], X_design_sim[!na_mask_from_output, ])), 
                 label = "X_whitened should differ from X_design_sim if global AR is non-zero (excluding initial NA rows)")
  } else {
    # If global AR coeffs are NULL or all zero, X_whitened should be X_design_sim
    expect_identical(X_w, X_design_sim, label = "X_whitened should be original X if global AR coeffs are null or zero.")
  }
})

test_that("ndx_ar2_whitening handles failed AR fits gracefully", {
  set.seed(456)
  n_voxels <- 2
  n_timepoints <- 100
  ar_order <- 2

  # Voxel 1: Normal AR data
  true_ar1_coeffs <- c(0.5, 0.1)
  ar_data1 <- .generate_ar_data(n_timepoints, true_ar1_coeffs, 1.0)
  
  # Voxel 2: Zero variance residuals (should cause AR fit to fail or produce NAs/zeros)
  zero_var_residuals <- rep(0, n_timepoints)

  Y_data_sim <- cbind(ar_data1$series, rnorm(n_timepoints)) # Y_data can be anything
  Y_residuals_for_AR_fit_sim <- cbind(ar_data1$series, zero_var_residuals)
  X_design_sim <- matrix(rnorm(n_timepoints * 2), n_timepoints, 2)

  whitening_results <- ndx_ar2_whitening(Y_data_sim, X_design_sim, Y_residuals_for_AR_fit_sim, order = ar_order)

  # Check AR_coeffs_voxelwise for the failed voxel
  expect_equal(whitening_results$AR_coeffs_voxelwise[2, ], c(0,0), 
               label = "AR coeffs for failed fit (zero var resid) should be c(0,0)")
  expect_true(is.na(whitening_results$var_innovations_voxelwise[2]), 
              label = "Innovation variance for failed fit should be NA")

  # Y_whitened for the failed voxel should be same as original Y_data for that voxel
  # (after accounting for NA rows if global_ar_on_design also resulted in non-filtering for X)
  # Current implementation means Y_whitened[,2] will be Y_data_sim[,2] as AR coeffs are 0
  na_mask_failed_fit <- whitening_results$na_mask
  expect_equal(whitening_results$Y_whitened[!na_mask_failed_fit, 2], Y_data_sim[!na_mask_failed_fit, 2], 
               label = "Y_whitened for failed AR fit voxel should be original data (excluding initial NA rows)")
  
  # Global AR coeffs should be based only on Voxel 1
  expect_true(!is.null(whitening_results$AR_coeffs_global))
  expect_equal(whitening_results$AR_coeffs_global, whitening_results$AR_coeffs_voxelwise[1,], tolerance = 1e-6,
               label = "Global AR coeffs should be based on the successful voxel fit")
})


test_that("ndx_ar2_whitening with global_ar_on_design = FALSE", {
  set.seed(789)
  n_voxels <- 1
  n_timepoints <- 50
  ar_order <- 1
  true_ar_coeffs <- c(0.7)
  ar_data <- .generate_ar_data(n_timepoints, true_ar_coeffs, 1.0)
  
  Y_data_sim <- matrix(ar_data$series, ncol = 1)
  Y_residuals_for_AR_fit_sim <- Y_data_sim
  X_design_sim <- matrix(rnorm(n_timepoints * 2), n_timepoints, 2)

  whitening_results <- ndx_ar2_whitening(Y_data_sim, X_design_sim, Y_residuals_for_AR_fit_sim, 
                                        order = ar_order, global_ar_on_design = FALSE)
  
  expect_true(is.list(whitening_results))
  expect_null(whitening_results$AR_coeffs_global, label = "AR_coeffs_global should be NULL when global_ar_on_design is FALSE")
  expect_identical(whitening_results$X_whitened, X_design_sim, 
                   label = "X_whitened should be identical to X_design_full when global_ar_on_design is FALSE")
  expect_true(!is.null(whitening_results$Y_whitened)) # Y_whitened should still be processed
})

test_that("Input validation for ndx_ar2_whitening", {
    Y_good <- matrix(rnorm(100*2), 100, 2)
    X_good <- matrix(rnorm(100*3), 100, 3)
    Y_res_good <- matrix(rnorm(100*2), 100, 2)

    expect_error(ndx_ar2_whitening(as.data.frame(Y_good), X_good, Y_res_good), "Y_data must be a numeric matrix.")
    expect_error(ndx_ar2_whitening(Y_good, as.data.frame(X_good), Y_res_good), "X_design_full must be a numeric matrix.")
    expect_error(ndx_ar2_whitening(Y_good, X_good, as.data.frame(Y_res_good)), "Y_residuals_for_AR_fit must be a numeric matrix.")

    Y_short_row <- matrix(rnorm(90*2), 90, 2)
    expect_error(ndx_ar2_whitening(Y_short_row, X_good, Y_res_good), "Y_data, X_design_full, and Y_residuals_for_AR_fit must have the same number of rows")
    
    Y_res_wrong_col <- matrix(rnorm(100*1), 100, 1)
    expect_error(ndx_ar2_whitening(Y_good, X_good, Y_res_wrong_col), "Y_data and Y_residuals_for_AR_fit must have the same number of columns")

    expect_error(ndx_ar2_whitening(Y_good, X_good, Y_res_good, order = 0.5), "order must be a single positive integer.")
    expect_error(ndx_ar2_whitening(Y_good, X_good, Y_res_good, order = 0), "order must be a single positive integer.")
    expect_error(ndx_ar2_whitening(Y_good, X_good, Y_res_good, order = c(1,2)), "order must be a single positive integer.")
    
    Y_too_short <- matrix(rnorm(2*2), 2, 2)
    X_too_short <- matrix(rnorm(2*3), 2, 3)
    Y_res_too_short <- matrix(rnorm(2*2), 2, 2)
    expect_error(ndx_ar2_whitening(Y_too_short, X_too_short, Y_res_too_short, order = 2L), 
                 regexp = "Number of timepoints \\(2\\) must be greater than AR order \\(2\\)\\.")
}) 

test_that("ndx_ar2_whitening effectively whitens AR(2) noise", {
  set.seed(2024)
  n_timepoints <- 400
  ar_coeffs <- c(0.8, -0.25)
  sim <- .generate_ar_data(n_timepoints, ar_coeffs, 1.0)
  Y <- matrix(sim$series, ncol = 1)
  X <- matrix(rnorm(n_timepoints * 2), n_timepoints, 2)
  res <- ndx_ar2_whitening(Y, X, Y, order = 2L, verbose = FALSE)
  Yw <- res$Y_whitened[!res$na_mask, 1]
  lb <- Box.test(Yw, lag = 5, type = "Ljung-Box")
  expect_gt(lb$p.value, 0.05, label = "Whitened series should pass Ljung-Box test")
  ar_after <- stats::ar.yw(Yw, aic = FALSE, order.max = 2)
  expect_true(all(abs(ar_after$ar) < 0.1), label = "AR coeffs after whitening near zero")
})

test_that("ndx_ar2_whitening leaves white noise unchanged", {
  set.seed(2025)
  n_timepoints <- 200
  Y <- matrix(rnorm(n_timepoints), ncol = 1)
  X <- matrix(rnorm(n_timepoints), ncol = 1)
  res <- ndx_ar2_whitening(Y, X, Y, order = 2L, verbose = FALSE)
  expect_true(all(abs(res$AR_coeffs_voxelwise) < 0.1), label = "AR coeffs of white noise near zero")
  expect_equal(res$Y_whitened[!res$na_mask, 1], Y[!res$na_mask, 1], tolerance = 0.05)
})
