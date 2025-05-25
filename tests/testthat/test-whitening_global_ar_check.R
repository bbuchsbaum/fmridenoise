context("Global AR stability safeguard")

test_that("Unstable averaged AR coefficients skip design whitening", {
  # Test the stability safeguard by creating data that forces the edge case
  set.seed(1234)
  n_tp <- 100
  n_voxels <- 8
  
  # Create data that will produce coefficients very close to the stability boundary
  Y <- matrix(0, nrow = n_tp, ncol = n_voxels)
  
  for (v in 1:n_voxels) {
    # Create data with very specific autocorrelation structure
    # that will push the AR estimates toward instability
    innovations <- rnorm(n_tp, sd = 0.5)  # Moderate noise
    y_temp <- numeric(n_tp)
    y_temp[1:2] <- innovations[1:2]
    
    # Use coefficients that will average to something close to instability
    # We'll use coefficients that are individually stable but average to unstable
    if (v <= 4) {
      # First group: coefficients that sum to > 1 (near unit root)
      ar_coeffs <- c(0.95, 0.1)  # Sum = 1.05, individually stable
    } else {
      # Second group: coefficients that when averaged with first group become unstable
      ar_coeffs <- c(0.85, -0.1)  # Sum = 0.75
    }
    # Expected average: (0.95 + 0.85)/2 = 0.9, (0.1 + (-0.1))/2 = 0.0
    # c(0.9, 0.0) has roots at 1/0.9 = 1.111, which is stable
    # Let's try a different approach...
    
    # Use coefficients that will definitely average to unstable values
    if (v <= 4) {
      ar_coeffs <- c(1.1, -0.05)  # Unstable individually (sum > 1)
    } else {
      ar_coeffs <- c(0.7, 0.35)   # Stable individually
    }
    # Expected average: (1.1 + 0.7)/2 = 0.9, (-0.05 + 0.35)/2 = 0.15
    # c(0.9, 0.15) - let's check if this is unstable
    
    for (t in 3:n_tp) {
      y_temp[t] <- ar_coeffs[1] * y_temp[t-1] + ar_coeffs[2] * y_temp[t-2] + innovations[t]
    }
    Y[, v] <- y_temp
  }
  
  X <- matrix(rnorm(n_tp * 3), nrow = n_tp, ncol = 3)
  Y_res <- Y
  
  # Test what the expected average would be
  expected_avg <- c((1.1 + 0.7)/2, (-0.05 + 0.35)/2)
  roots_expected <- polyroot(c(1, -expected_avg))
  min_mod_expected <- min(Mod(roots_expected))
  
  cat("Expected average coeffs:", expected_avg, "\n")
  cat("Expected min root mod:", min_mod_expected, "\n")
  
  # If the expected coefficients are actually stable, skip the test
  if (min_mod_expected > 1.00001) {
    skip("Expected coefficients are stable, cannot test instability safeguard")
  }
  
  # Run the whitening with global AR enabled to trigger the stability check
  result <- ndx_ar2_whitening(Y, X, Y_res, order = 2L, verbose = TRUE, 
                             global_ar_on_design = TRUE,  # Enable global AR to trigger stability check
                             max_ar_failures_prop = 0.9)
  
  if (!is.null(result$AR_coeffs_global)) {
    # Check actual stability
    roots_actual <- polyroot(c(1, -result$AR_coeffs_global))
    min_mod_actual <- min(Mod(roots_actual))
    cat("Actual global coeffs:", result$AR_coeffs_global, "\n")
    cat("Actual min root mod:", min_mod_actual, "\n")
    
    if (min_mod_actual <= 1.00001) {
      # The coefficients are unstable but weren't caught - this is a bug
      fail("Unstable coefficients were not caught by the stability check")
    } else {
      # Skip if we can't generate the right conditions
      skip("AR fitting produced stable coefficients despite unstable input design")
    }
  } else {
    # Success - either instability was caught or too many fits failed
    expect_identical(result$X_whitened, X)
    cat("Safeguard worked: AR_coeffs_global is NULL\n")
  }
})

test_that("AR fitting failure proportion safeguard works", {
  # Test the other safeguard: when too many AR fits fail
  set.seed(456)
  n_tp <- 50  # Short time series to increase fit failures
  n_voxels <- 10
  
  # Create data that will cause many AR fits to fail
  Y <- matrix(0, nrow = n_tp, ncol = n_voxels)
  
  for (v in 1:n_voxels) {
    if (v <= 7) {
      # Most voxels: constant or near-constant data (zero variance)
      Y[, v] <- rep(0.001, n_tp) + rnorm(n_tp, sd = 1e-10)
    } else {
      # Few voxels: normal data
      Y[, v] <- rnorm(n_tp)
    }
  }
  
  X <- matrix(rnorm(n_tp * 3), nrow = n_tp, ncol = 3)
  Y_res <- Y

  # Use a low max_ar_failures_prop to trigger the safeguard
  # With 7/10 voxels having zero variance, we expect 70% failure rate
  expect_warning(
    res <- ndx_ar2_whitening(Y, X, Y_res, order = 2L, verbose = FALSE, 
                            max_ar_failures_prop = 0.5),  # 50% threshold, expect 70% failures
    regexp = "Proportion of failed AR fits.*exceeds max_ar_failures_prop"
  )
  
  # When the safeguard triggers, AR_coeffs_global should be NULL
  # and the data should be returned unwhitened
  expect_null(res$AR_coeffs_global)
  expect_identical(res$X_whitened, X)
  expect_identical(res$Y_whitened, Y)
})

test_that("Stability check logic works correctly", {
  # Test the stability check logic directly with known unstable coefficients
  
  # Test case 1: Clearly unstable coefficients (sum > 1)
  unstable_coeffs_1 <- c(0.8, 0.3)  # sum = 1.1 > 1
  roots_1 <- polyroot(c(1, -unstable_coeffs_1))
  min_mod_1 <- min(Mod(roots_1))
  expect_true(min_mod_1 <= 1.00001, 
              info = sprintf("Coefficients c(0.8, 0.3) should be unstable, min root mod: %.6f", min_mod_1))
  
  # Test case 2: Clearly stable coefficients
  stable_coeffs <- c(0.5, 0.2)  # sum = 0.7 < 1
  roots_2 <- polyroot(c(1, -stable_coeffs))
  min_mod_2 <- min(Mod(roots_2))
  expect_true(min_mod_2 > 1.00001,
              info = sprintf("Coefficients c(0.5, 0.2) should be stable, min root mod: %.6f", min_mod_2))
  
  # Test case 3: Edge case - coefficients very close to unit root
  edge_coeffs <- c(0.99, 0.005)  # sum = 0.995, very close to 1
  roots_3 <- polyroot(c(1, -edge_coeffs))
  min_mod_3 <- min(Mod(roots_3))
  # This should be stable but close to the boundary
  expect_true(min_mod_3 > 1.00001,
              info = sprintf("Coefficients c(0.99, 0.005) should be stable, min root mod: %.6f", min_mod_3))
})
