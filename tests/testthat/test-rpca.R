context("ndx_rpca_temporal_components_multirun - Functionality")

# Helper to generate some Y_residuals_cat and run_idx
.generate_rpca_test_data <- function(T_run = 100, V = 50, N_runs = 1) {
  total_T <- T_run * N_runs
  Y_res_cat <- matrix(rnorm(total_T * V), total_T, V) * 0.01
  run_idx_vec <- rep(1:N_runs, each = T_run)
  
  if (V >= 2) { 
      true_V_pattern <- matrix(rnorm(V*2), V, 2) 
      for (r in 1:N_runs) {
        run_rows <- which(run_idx_vec == r)
        # Stronger low-rank signal
        C_r_signal_run1 <- sin((1:T_run)/8 + r/2) * 5 
        C_r_signal_run2 <- cos((1:T_run)/15 - r/3) * 4 
        Y_res_cat[run_rows, ] <- Y_res_cat[run_rows, ] + 
                                   (cbind(C_r_signal_run1, C_r_signal_run2) %*% t(true_V_pattern)) * 5.0 # Increased signal multiplier
      }
  }
  return(list(Y_residuals_cat = Y_res_cat, run_idx = run_idx_vec, T_run=T_run, V=V, N_runs=N_runs, total_T=total_T))
}

test_that("ndx_rpca_temporal_components_multirun runs (single run, concat_svd)", {
  test_data <- .generate_rpca_test_data(N_runs = 1)
  k_target <- 3
  
  user_opts <- list(
    k_per_run_target = 4, 
    rpca_lambda_auto = TRUE, # Revert to auto lambda
    # rpca_lambda_fixed = 0.05, # Try a fixed, possibly smaller lambda
    rpca_mu = 0.3, # Try an explicit mu for test data
    rpca_merge_strategy = "concat_svd",
    rpca_term_delta = 1e-6, # Explicitly use function's default
    rpca_max_iter = 2000    # Explicitly use function's default
  )
  
  components <- NULL
  expect_no_error({
    components <- ndx_rpca_temporal_components_multirun(
      Y_residuals_cat = test_data$Y_residuals_cat,
      run_idx = test_data$run_idx,
      k_global_target = k_target,
      user_options = user_opts
    )
  })
  
  expect_true(!is.null(components), "Components should not be NULL for valid inputs (single run)")
  if (!is.null(components)) {
    expect_true(is.matrix(components), "Output should be a matrix (single run)")
    expect_equal(nrow(components), test_data$total_T, info = "Output rows should match total timepoints (single run)")
    expect_true(ncol(components) <= k_target, info = sprintf("Output columns (%d) should be <= k_target (%d) (single run)", ncol(components), k_target))
    expect_true(ncol(components) > 0, info = "Output should have >0 components if k_target > 0 (single run)")
  }
})

test_that("ndx_rpca_temporal_components_multirun runs (multi-run, concat_svd)", {
  test_data <- .generate_rpca_test_data(N_runs = 2)
  k_target <- 3
  
  user_opts <- list(
    k_per_run_target = 4,
    rpca_lambda_auto = TRUE, # Revert to auto lambda
    # rpca_lambda_fixed = 0.05, # Try a fixed, possibly smaller lambda
    rpca_mu = 0.3, # Try an explicit mu for test data
    rpca_merge_strategy = "concat_svd",
    rpca_term_delta = 1e-6, # Explicitly use function's default
    rpca_max_iter = 2000    # Explicitly use function's default
  )
  
  components <- NULL
  expect_no_error({
    components <- ndx_rpca_temporal_components_multirun(
      Y_residuals_cat = test_data$Y_residuals_cat,
      run_idx = test_data$run_idx,
      k_global_target = k_target,
      user_options = user_opts
    )
  })
  
  expect_true(!is.null(components), "Components should not be NULL for valid inputs (multi-run)")
  if (!is.null(components)) {
    expect_true(is.matrix(components), "Output should be a matrix (multi-run)")
    expect_equal(nrow(components), test_data$total_T, info = "Output rows should match total timepoints (multi-run)")
    expect_true(ncol(components) <= k_target, info = sprintf("Output columns (%d) should be <= k_target (%d) (multi-run)", ncol(components), k_target))
    expect_true(ncol(components) > 0, info = "Output should have >0 components if k_target > 0 (multi-run)")
  }
}) 