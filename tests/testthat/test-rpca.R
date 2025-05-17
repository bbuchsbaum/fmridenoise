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
    expect_true(is.list(components))
    expect_true(is.matrix(components$C_components))
    expect_equal(nrow(components$C_components), test_data$total_T)
    expect_true(ncol(components$C_components) <= k_target)
    expect_length(components$spike_TR_mask, test_data$total_T)
    expect_type(components$spike_TR_mask, "logical")
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
    expect_true(is.list(components))
    expect_true(is.matrix(components$C_components))
    expect_equal(nrow(components$C_components), test_data$total_T)
    expect_true(ncol(components$C_components) <= k_target)
    expect_length(components$spike_TR_mask, test_data$total_T)
    expect_type(components$spike_TR_mask, "logical")
  }
})

test_that("ndx_rpca_temporal_components_multirun runs with iterative merge strategy", {
  test_data <- .generate_rpca_test_data(N_runs = 3)
  k_target <- 2

  user_opts <- list(
    k_per_run_target = 3,
    rpca_lambda_auto = TRUE,
    rpca_mu = 0.3,
    rpca_merge_strategy = "iterative",
    rpca_term_delta = 1e-6,
    rpca_max_iter = 2000
  )

  components_iter <- NULL
  expect_no_error({
    components_iter <- ndx_rpca_temporal_components_multirun(
      Y_residuals_cat = test_data$Y_residuals_cat,
      run_idx = test_data$run_idx,
      k_global_target = k_target,
      user_options = user_opts
    )
  })

  expect_true(!is.null(components_iter), "Components should not be NULL for iterative merge")
  if (!is.null(components_iter)) {
    expect_true(is.list(components_iter))
    expect_true(is.matrix(components_iter$C_components))
    expect_equal(nrow(components_iter$C_components), test_data$total_T)
    expect_true(ncol(components_iter$C_components) <= k_target)
    expect_length(components_iter$spike_TR_mask, test_data$total_T)
    expect_type(components_iter$spike_TR_mask, "logical")

    # The C_r components (components_iter) are not guaranteed to be orthonormal.
    # V_global (the voxel-space basis) is orthonormal, but C_r = E_r %*% V_global.
    # Removing orthonormality checks on C_r directly.
    # if (ncol(components_iter) > 1) { # only check if more than 1 component
    #   gram <- crossprod(components_iter)
    #   off_diag_max <- max(abs(gram[upper.tri(gram)]))
    #   diag_dev_max <- max(abs(diag(gram) - 1))
    #   expect_lt(off_diag_max, 1e-6, 
    #             label = "Iterative merge components should be nearly orthogonal")
    #   expect_lt(diag_dev_max, 1e-6, 
    #             label = "Iterative merge components should have unit norm")
    # }
  }
})
