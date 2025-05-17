context("ndx_initial_glm - Basic Functionality")

# Mock data for testing ndx_initial_glm
# Similar setup to what might be in ndx_run_sprint1 for this step

test_that("ndx_initial_glm runs with minimal valid inputs and returns correct structure", {
  set.seed(123) # For reproducibility
  TR_example <- 2.0
  n_time_per_run <- 10
  n_runs <- 1
  n_voxels <- 3
  total_timepoints <- n_time_per_run * n_runs

  Y_fmri_example <- matrix(rnorm(total_timepoints * n_voxels), nrow = total_timepoints, ncol = n_voxels)
  run_idx_example <- rep(1:n_runs, each = n_time_per_run)
  motion_params_example <- matrix(rnorm(total_timepoints * 6), nrow = total_timepoints, ncol = 6)
  colnames(motion_params_example) <- paste0("mot", 1:6) # Add colnames to motion_params

  events_example <- data.frame(
    onsets = as.numeric(c(2, 5)),
    durations = as.numeric(c(1, 1)),
    condition = factor(c("TaskA", "TaskB")),
    blockids = as.integer(c(1, 1))
  )
  
  # It's good practice to wrap the call in a tryCatch or expect_no_error
  # if the main point is just that it runs and has structure.
  # For now, we just expect it to run.
  output <- NULL
  expect_no_error({
      output <- ndx_initial_glm(Y_fmri_example, events_example, 
                                   motion_params_example, run_idx_example, TR_example)
  })
  
  # Basic checks
  expect_true(is.list(output), "Output should be a list")
  expect_named(output, c("Y_residuals_current", "VAR_BASELINE_FOR_DES"))
  expect_true(is.matrix(output$Y_residuals_current))
  expect_equal(dim(output$Y_residuals_current), dim(Y_fmri_example))
  expect_true(is.numeric(output$VAR_BASELINE_FOR_DES) && length(output$VAR_BASELINE_FOR_DES) == 1)
  expect_false(any(is.na(output$Y_residuals_current))) # OLS residuals should be complete if inputs are
  expect_true(is.finite(output$VAR_BASELINE_FOR_DES))
})

test_that("ndx_initial_glm errors if Y_fmri is not a matrix", {
  set.seed(123)
  TR_example <- 2.0
  n_time_per_run <- 10
  n_runs <- 1
  n_voxels <- 3
  total_timepoints <- n_time_per_run * n_runs

  Y_fmri_example <- matrix(rnorm(total_timepoints * n_voxels), nrow = total_timepoints, ncol = n_voxels)
  run_idx_example <- rep(1:n_runs, each = n_time_per_run)
  motion_params_example <- matrix(rnorm(total_timepoints * 6), nrow = total_timepoints, ncol = 6)
  colnames(motion_params_example) <- paste0("mot", 1:6)

  events_example <- data.frame(
    onsets = as.numeric(c(2, 5)),
    durations = as.numeric(c(1, 1)),
    condition = factor(c("TaskA", "TaskB")),
    blockids = as.integer(c(1, 1))
  )

  Y_fmri_invalid <- as.data.frame(Y_fmri_example)
  
  expect_error(
    ndx_initial_glm(Y_fmri_invalid, events_example,
                    motion_params_example, run_idx_example, TR_example),
    "Y_fmri must be a matrix"
  )
})

test_that("ndx_initial_glm errors if events data frame is missing required columns", {
  set.seed(123)
  TR_example <- 2.0
  n_time_per_run <- 10
  n_runs <- 1
  n_voxels <- 3
  total_timepoints <- n_time_per_run * n_runs

  Y_fmri_example <- matrix(rnorm(total_timepoints * n_voxels), nrow = total_timepoints, ncol = n_voxels)
  run_idx_example <- rep(1, total_timepoints)
  motion_params_example <- matrix(rnorm(total_timepoints * 6), nrow = total_timepoints, ncol = 6)
  colnames(motion_params_example) <- paste0("mot", 1:6)

  events_example <- data.frame(
    onsets = as.numeric(c(2, 5)),
    durations = as.numeric(c(1, 1)),
    condition = factor(c("TaskA", "TaskB")),
    blockids = as.integer(c(1, 1))
  )

  events_invalid <- events_example[, -which(names(events_example) == "onsets")]
  expect_error(
    ndx_initial_glm(Y_fmri_example, events_invalid,
                    motion_params_example, run_idx_example, TR_example),
    "events data frame must contain columns: onsets, durations, condition, blockids"
  )
  
  events_invalid_cond <- events_example
  events_invalid_cond$condition <- NULL
  expect_error(
    ndx_initial_glm(Y_fmri_example, events_invalid_cond,
                    motion_params_example, run_idx_example, TR_example),
    "events data frame must contain columns: onsets, durations, condition, blockids"
  )
})

test_that("ndx_initial_glm handles multiple runs correctly", {
  set.seed(456)
  TR_multirun <- 2.5
  n_time_per_run <- 10
  n_runs <- 2
  n_voxels <- 3
  total_timepoints <- n_time_per_run * n_runs

  Y_fmri_multirun <- matrix(rnorm(total_timepoints * n_voxels), nrow = total_timepoints, ncol = n_voxels)
  run_idx_multirun <- rep(1:n_runs, each = n_time_per_run)
  motion_params_multirun <- matrix(rnorm(total_timepoints * 6), nrow = total_timepoints, ncol = 6)
  colnames(motion_params_multirun) <- paste0("mot", 1:6)

  events_multirun <- data.frame(
    onsets = as.numeric(c(2, 5, 2, 5)), # Events in run 1, then run 2
    durations = as.numeric(c(1, 1, 1, 1)),
    condition = factor(rep(c("TaskA", "TaskB"), 2)),
    blockids = as.integer(c(1, 1, 2, 2)) # Correct blockids for each run
  )
  
  output_multirun <- NULL
  expect_no_error({
    output_multirun <- ndx_initial_glm(Y_fmri_multirun, events_multirun, 
                                         motion_params_multirun, run_idx_multirun, TR_multirun)
  })
  
  # Check structure (same as single run)
  expect_true(is.list(output_multirun), "Output should be a list")
  expect_named(output_multirun, c("Y_residuals_current", "VAR_BASELINE_FOR_DES"))
  expect_true(is.matrix(output_multirun$Y_residuals_current))
  expect_equal(dim(output_multirun$Y_residuals_current), dim(Y_fmri_multirun))
  expect_true(is.numeric(output_multirun$VAR_BASELINE_FOR_DES) && length(output_multirun$VAR_BASELINE_FOR_DES) == 1)
  expect_false(any(is.na(output_multirun$Y_residuals_current)))
  expect_true(is.finite(output_multirun$VAR_BASELINE_FOR_DES))
  
  # Further checks could compare to manual calculation or known results if simple enough.
  # For now, primarily checking it runs without error and structure is okay.
})

# Add more tests as needed: e.g., specific scenarios for VAR_BASELINE_FOR_DES, 
# edge cases with motion parameters, different polynomial degrees if that becomes an option, etc. 