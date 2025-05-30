context("ndx_build_design_matrix - Design Matrix Construction")

# --- Helper function to create mock estimated_hrfs tibble ---
create_mock_hrfs <- function(conditions = c("TaskA", "TaskB"), taps_per_hrf = 6, TR = 2.0) {
  if (length(conditions) == 0) return(NULL)
  hrf_list <- lapply(conditions, function(cond) {
    list(
      condition = cond,
      hrf_estimate = list(stats::rnorm(taps_per_hrf)), # list-column
      taps = list(1:taps_per_hrf) # list-column, or could be just num_taps
    )
  })
  # Convert list of lists to a tibble more carefully
  # Each element of hrf_list is a row
  df_hrfs <- do.call(rbind, lapply(hrf_list, function(row_list) data.frame(condition=row_list$condition)))
  df_hrfs$hrf_estimate <- lapply(hrf_list, function(row_list) row_list$hrf_estimate[[1]]) # unlist the inner list for hrf_estimate
  df_hrfs$taps <- lapply(hrf_list, function(row_list) row_list$taps[[1]]) # unlist the inner list for taps
  return(tibble::as_tibble(df_hrfs))
}

# --- Basic Test Data Setup ---
TR_test_dm <- 2.0
n_time_per_run_dm <- 30 
n_runs_dm <- 2
total_timepoints_dm <- n_time_per_run_dm * n_runs_dm
run_idx_dm <- rep(1:n_runs_dm, each = n_time_per_run_dm)

events_dm <- data.frame(
  onsets = c(5, 15, 5, 20) * TR_test_dm, # in seconds
  durations = rep(2 * TR_test_dm, 4),
  condition = factor(c("TaskA", "TaskB", "TaskA", "TaskB")),
  blockids = c(1, 1, 2, 2) # ensure blockids align with runs
)

motion_params_dm <- matrix(stats::rnorm(total_timepoints_dm * 3), ncol = 3)
colnames(motion_params_dm) <- paste0("mot", 1:3)

rpca_comps_dm <- matrix(stats::rnorm(total_timepoints_dm * 2), ncol = 2)
# colnames for rpca will be auto-generated by ndx_build_design_matrix

spectral_sines_dm <- matrix(stats::rnorm(total_timepoints_dm * 4), ncol = 4)
colnames(spectral_sines_dm) <- paste0(rep(c("s1", "s2"), each=2), c("_sin", "_cos"))

estimated_hrfs_dm <- create_mock_hrfs(conditions = c("TaskA", "TaskB"), taps_per_hrf = 8, TR = TR_test_dm)

# --- Test Cases ---
test_that("ndx_build_design_matrix runs with all components present (multi-run)", {
  X_full <- NULL
  expect_no_error({
    X_full <- ndx_build_design_matrix(
      estimated_hrfs = estimated_hrfs_dm,
      events = events_dm,
      motion_params = motion_params_dm,
      rpca_components = rpca_comps_dm,
      spectral_sines = spectral_sines_dm,
      run_idx = run_idx_dm,
      TR = TR_test_dm,
      poly_degree = 1, # poly0, poly1 for degree=1 (3 columns total: 1 global + 1 per run)
      verbose = FALSE
    )
  })
  
  expect_true(is.matrix(X_full))
  expect_equal(nrow(X_full), total_timepoints_dm)
  
  # Expected columns: 
  # TaskA (1), TaskB (1) = 2
  # Motion (3) = 3
  # RPCA (2) = 2 
  # Spectral (4) = 4
  # Poly (poly0, poly1_block_1, poly1_block_2 for degree=1) = 3
  # Run Intercepts (for run2, since poly0 is overall intercept and 2 runs) = 1 
  # Total = 2 + 3 + 2 + 4 + 3 + 1 = 15
  expect_equal(ncol(X_full), 15)
  
  expected_colnames_structure <- c(
    "task_TaskA", "task_TaskB", 
    paste0("mot", 1:3), 
    paste0("rpca_comp_", 1:2), 
    colnames(spectral_sines_dm), 
    "poly0", "poly1", "poly2",
    "run_intercept_2"
  )
  expect_equal(sort(colnames(X_full)), sort(expected_colnames_structure), 
               info = "Column names for X_full do not match expected structure (set comparison).")
})

test_that("ndx_build_design_matrix handles NULL/empty optional components", {
  X_no_nuisance <- NULL
  expect_no_error({
    X_no_nuisance <- ndx_build_design_matrix(
      estimated_hrfs = estimated_hrfs_dm,
      events = events_dm,
      motion_params = NULL,
      rpca_components = NULL,
      spectral_sines = NULL,
      run_idx = run_idx_dm,
      TR = TR_test_dm,
      poly_degree = 0, # only poly0 (overall intercept)
      verbose = FALSE
    )
  })
  expect_true(is.matrix(X_no_nuisance))
  expect_equal(nrow(X_no_nuisance), total_timepoints_dm)
  expect_equal(ncol(X_no_nuisance), 4)
  expected_colnames_no_nuisance <- c("task_TaskA", "task_TaskB", "poly0", "run_intercept_2")
  expect_equal(sort(colnames(X_no_nuisance)), sort(expected_colnames_no_nuisance),
               info = "Colnames for X_no_nuisance do not match (set comparison).")
  
  # Test with 0-column matrices for optional components
  X_zero_col_nuisance <- NULL
  expect_no_error({
    X_zero_col_nuisance <- ndx_build_design_matrix(
      estimated_hrfs = estimated_hrfs_dm,
      events = events_dm,
      motion_params = matrix(numeric(0), nrow=total_timepoints_dm, ncol=0),
      rpca_components = matrix(numeric(0), nrow=total_timepoints_dm, ncol=0),
      spectral_sines = matrix(numeric(0), nrow=total_timepoints_dm, ncol=0),
      run_idx = run_idx_dm,
      TR = TR_test_dm,
      poly_degree = -1, # No polynomials, only run intercepts if multi-run
      verbose = FALSE
    )
  })
  expect_true(is.matrix(X_zero_col_nuisance))
  expect_equal(ncol(X_zero_col_nuisance), 4)
  expected_colnames_zero_col_nuisance <- c("task_TaskA", "task_TaskB", "run_intercept_1", "run_intercept_2")
  expect_equal(sort(colnames(X_zero_col_nuisance)), sort(expected_colnames_zero_col_nuisance),
               info = "Colnames for X_zero_col_nuisance do not match (set comparison).")
})

test_that("ndx_build_design_matrix handles single run correctly", {
  run_idx_single_dm <- rep(1, n_time_per_run_dm)
  events_single_dm <- events_dm[events_dm$blockids == 1,]
  motion_single_dm <- motion_params_dm[1:n_time_per_run_dm, , drop=FALSE]
  rpca_single_dm <- rpca_comps_dm[1:n_time_per_run_dm, , drop=FALSE]
  spectral_single_dm <- spectral_sines_dm[1:n_time_per_run_dm, , drop=FALSE]
  
  X_single_run <- NULL
  expect_no_error({
    X_single_run <- ndx_build_design_matrix(
      estimated_hrfs = estimated_hrfs_dm, 
      events = events_single_dm,
      motion_params = motion_single_dm,
      rpca_components = rpca_single_dm,
      spectral_sines = spectral_single_dm,
      run_idx = run_idx_single_dm,
      TR = TR_test_dm,
      poly_degree = 1, # poly0, poly1_block_1 for degree=1 (2 columns for single run)
      verbose = FALSE
    )
  })
  expect_true(is.matrix(X_single_run))
  expect_equal(nrow(X_single_run), n_time_per_run_dm)
  # Task (2) + Motion (3) + RPCA (2) + Spectral (4) + Poly (2: poly0, poly1_block_1 for degree=1 on single run) = 13. 
  # No run-specific intercept as only 1 run and poly0 exists.
  expect_equal(ncol(X_single_run), 13)
  expect_true(all(c("poly0", "poly1") %in% colnames(X_single_run)))
  
  # Single run, no polynomials -> should add an overall intercept
  X_single_run_no_poly <- NULL
  expect_no_error({
    X_single_run_no_poly <- ndx_build_design_matrix(
      estimated_hrfs = estimated_hrfs_dm,
      events = events_single_dm,
      motion_params = motion_single_dm,
      rpca_components = NULL,
      spectral_sines = NULL,
      run_idx = run_idx_single_dm,
      TR = TR_test_dm,
      poly_degree = -1, # No polys
      verbose = FALSE
    )
  })
  expect_true(is.matrix(X_single_run_no_poly))
  # Task (2) + Motion (3) + Intercept (1) = 6
  expect_equal(ncol(X_single_run_no_poly), 6)
  expect_true("intercept" %in% colnames(X_single_run_no_poly))
})

test_that("ndx_build_design_matrix handles no task HRFs", {
  X_no_task <- NULL
  expect_no_error({
    X_no_task <- ndx_build_design_matrix(
      estimated_hrfs = NULL,
      events = events_dm,
      motion_params = motion_params_dm,
      rpca_components = rpca_comps_dm,
      spectral_sines = spectral_sines_dm,
      run_idx = run_idx_dm,
      TR = TR_test_dm,
      poly_degree = 1, # poly0, poly1_block_1, poly1_block_2
      verbose = FALSE
    )
  })
  expect_true(is.matrix(X_no_task))
  # Expected: Motion (3) + RPCA (2) + Spectral (4) + Poly (3 for degree=1) + Run Intercept (1) = 13
  expect_equal(ncol(X_no_task), 13)
  expect_false(any(startsWith(colnames(X_no_task), "task_")))
  
  # Empty tibble for HRFs
  empty_hrfs <- tibble::tibble(condition=character(), hrf_estimate=list(), taps=list())
  expect_no_error({
    X_no_task_empty_hrf <- ndx_build_design_matrix(
      estimated_hrfs = empty_hrfs,
      events = events_dm,
      motion_params = motion_params_dm,
      rpca_components = rpca_comps_dm,
      spectral_sines = spectral_sines_dm,
      run_idx = run_idx_dm,
      TR = TR_test_dm,
      poly_degree = 1, # poly0, poly1_block_1, poly1_block_2
      verbose = FALSE
    )
  })
  expect_equal(ncol(X_no_task_empty_hrf), 13)
})

test_that("ndx_build_design_matrix generates correct polynomial degrees and intercepts", {
  # poly_degree = 0 (intercept only from poly)
  X_poly0 <- ndx_build_design_matrix(estimated_hrfs=NULL, events=events_dm, motion_params=NULL, rpca_components=NULL, spectral_sines=NULL, 
                                   run_idx=run_idx_dm, TR=TR_test_dm, poly_degree=0, verbose=FALSE)
  # poly0 (1) + run_intercept_2 (1) = 2 columns
  expect_equal(ncol(X_poly0), 2)
  expect_true(all(c("poly0", "run_intercept_2") %in% colnames(X_poly0)))
  
  # poly_degree = -1 (no poly, run intercepts only if multi-run)
  X_no_poly_multirun <- ndx_build_design_matrix(estimated_hrfs=NULL, events=events_dm, motion_params=NULL, rpca_components=NULL, spectral_sines=NULL, 
                                             run_idx=run_idx_dm, TR=TR_test_dm, poly_degree= -1, verbose=FALSE)
  # run_intercept_1, run_intercept_2 = 2 columns
  expect_equal(ncol(X_no_poly_multirun), 2)
  expect_true(all(c("run_intercept_1", "run_intercept_2") %in% colnames(X_no_poly_multirun)))
  
  # poly_degree = -1 (no poly, single run -> should add overall intercept)
  run_idx_single_dm_for_poly <- rep(1, n_time_per_run_dm)
  X_no_poly_singlerun <- ndx_build_design_matrix(estimated_hrfs=NULL, events=events_dm[events_dm$blockids==1,], motion_params=NULL, rpca_components=NULL, spectral_sines=NULL, 
                                                run_idx=run_idx_single_dm_for_poly, TR=TR_test_dm, poly_degree= -1, verbose=FALSE)
  # intercept (1) = 1 column
  expect_equal(ncol(X_no_poly_singlerun), 1)
  expect_true("intercept" %in% colnames(X_no_poly_singlerun))
  
  # Test higher polynomial degrees
  X_poly2 <- ndx_build_design_matrix(estimated_hrfs=NULL, events=events_dm, motion_params=NULL, rpca_components=NULL, spectral_sines=NULL, 
                                   run_idx=run_idx_dm, TR=TR_test_dm, poly_degree=2, verbose=FALSE)
  # poly0, poly1_block_1, poly2_block_1, poly1_block_2, poly2_block_2 (5) + run_intercept_2 (1) = 6 columns
  expect_equal(ncol(X_poly2), 6)
  expect_true(all(c("poly0", "poly1", "poly2", "poly3", "poly4", "run_intercept_2") %in% colnames(X_poly2)))
})

test_that("ndx_build_design_matrix errors on row mismatches for inputs", {
  bad_motion <- motion_params_dm[1:(total_timepoints_dm-1), , drop=FALSE]
  expect_error(ndx_build_design_matrix(estimated_hrfs=NULL, events=events_dm, motion_params=bad_motion, rpca_components=NULL, spectral_sines=NULL,
                                      run_idx=run_idx_dm, TR=TR_test_dm, poly_degree=0, verbose=FALSE),
               "Row mismatch: motion_params")
  
  bad_rpca <- rpca_comps_dm[1:(total_timepoints_dm-2), , drop=FALSE]
  expect_error(ndx_build_design_matrix(estimated_hrfs=NULL, events=events_dm, motion_params=NULL, rpca_components=bad_rpca, spectral_sines=NULL,
                                      run_idx=run_idx_dm, TR=TR_test_dm, poly_degree=0, verbose=FALSE),
               "Row mismatch: rpca_components")
  
  bad_spectral <- spectral_sines_dm[1:(total_timepoints_dm-3), , drop=FALSE]
  expect_error(ndx_build_design_matrix(estimated_hrfs=NULL, events=events_dm, motion_params=NULL, rpca_components=NULL, spectral_sines=bad_spectral,
                                      run_idx=run_idx_dm, TR=TR_test_dm, poly_degree=0, verbose=FALSE),
               "Row mismatch: spectral_sines")
})

test_that("ndx_build_design_matrix errors on invalid TR", {
  expect_error(ndx_build_design_matrix(estimated_hrfs=NULL, events=events_dm, motion_params=NULL, rpca_components=NULL, spectral_sines=NULL,
                                      run_idx=run_idx_dm, TR=c(2, 3), poly_degree=0, verbose=FALSE),
               "single positive number")
  expect_error(ndx_build_design_matrix(estimated_hrfs=NULL, events=events_dm, motion_params=NULL, rpca_components=NULL, spectral_sines=NULL,
                                      run_idx=run_idx_dm, TR=-1, poly_degree=0, verbose=FALSE),
               "single positive number")
  expect_error(ndx_build_design_matrix(estimated_hrfs=NULL, events=events_dm, motion_params=NULL, rpca_components=NULL, spectral_sines=NULL,
                                      run_idx=run_idx_dm, TR=0, poly_degree=0, verbose=FALSE),
               "single positive number")
  expect_error(ndx_build_design_matrix(estimated_hrfs=NULL, events=events_dm, motion_params=NULL, rpca_components=NULL, spectral_sines=NULL,
                                      run_idx=run_idx_dm, TR=NA, poly_degree=0, verbose=FALSE),
               "single positive number")
})

test_that("ndx_build_design_matrix errors on invalid run_idx", {
  # Test with run_idx implying 0 timepoints
  expect_error(ndx_build_design_matrix(estimated_hrfs=NULL, events=events_dm, motion_params=NULL, rpca_components=NULL, spectral_sines=NULL, 
                                     run_idx=integer(0), TR=2, poly_degree=0, verbose=FALSE),
               "run_idx implies one or more runs have zero or negative length")
  
  # Test with NA values in run_idx
  bad_run_idx <- run_idx_dm
  bad_run_idx[5] <- NA
  expect_error(
    ndx_build_design_matrix(
      estimated_hrfs = estimated_hrfs_dm,
      events = events_dm,
      motion_params = motion_params_dm,
      rpca_components = rpca_comps_dm,
      spectral_sines = spectral_sines_dm,
      run_idx = bad_run_idx,
      TR = TR_test_dm,
      poly_degree = 1,
      verbose = FALSE
    ),
    "run_idx contains NA"
  )
  
  # Test with non-integer values in run_idx
  bad_run_idx_float <- run_idx_dm
  bad_run_idx_float[5] <- 1.5
  expect_error(
    ndx_build_design_matrix(
      estimated_hrfs = estimated_hrfs_dm,
      events = events_dm,
      motion_params = motion_params_dm,
      rpca_components = rpca_comps_dm,
      spectral_sines = spectral_sines_dm,
      run_idx = bad_run_idx_float,
      TR = TR_test_dm,
      poly_degree = 1,
      verbose = FALSE
    ),
    "run_idx must contain integer values only"
  )
  
  # Test with infinite values in run_idx
  bad_run_idx_inf <- run_idx_dm
  bad_run_idx_inf[5] <- Inf
  expect_error(
    ndx_build_design_matrix(
      estimated_hrfs = estimated_hrfs_dm,
      events = events_dm,
      motion_params = motion_params_dm,
      rpca_components = rpca_comps_dm,
      spectral_sines = spectral_sines_dm,
      run_idx = bad_run_idx_inf,
      TR = TR_test_dm,
      poly_degree = 1,
      verbose = FALSE
    ),
    "run_idx contains NA or non-finite values"
  )
})

# Test specific logic for FIR basis generation if complex conditions arise
test_that("FIR basis generation in ndx_build_design_matrix handles edge cases", {
  run_idx_single_dm <- rep(1, n_time_per_run_dm)
  events_single_dm <- events_dm[events_dm$blockids == 1,]
  
  # Case 1: HRF estimate has zero length
  hrf_zero_len <- create_mock_hrfs(conditions = c("TaskA"), taps_per_hrf = 0, TR = TR_test_dm)
  # This should be caught by: if (is.null(hrf_coeffs) || length(hrf_coeffs) == 0)
  X_fir_zero <- ndx_build_design_matrix(estimated_hrfs = hrf_zero_len, events = events_single_dm, motion_params = NULL, rpca_components = NULL, 
                                      spectral_sines = NULL, run_idx = run_idx_single_dm, TR = TR_test_dm, poly_degree = -1, verbose = FALSE)
  expect_false(any(startsWith(colnames(X_fir_zero), "task_")))
  expect_equal(ncol(X_fir_zero), 1) # Should be just the intercept

  # Case 2: No events for a condition listed in HRF table
  hrf_taskC <- create_mock_hrfs(conditions = c("TaskC"), taps_per_hrf = 6, TR = TR_test_dm)
  X_taskC_no_events <- ndx_build_design_matrix(estimated_hrfs = hrf_taskC, events = events_single_dm, motion_params = NULL, rpca_components = NULL,
                                             spectral_sines = NULL, run_idx = run_idx_single_dm, TR = TR_test_dm, poly_degree = -1, verbose = FALSE)
  expect_false(any(startsWith(colnames(X_taskC_no_events), "task_")))
  expect_equal(ncol(X_taskC_no_events), 1) # Intercept
  
  # Case 3: HRF estimate is NULL for a condition
  hrf_null_estimate <- tibble::tibble(
    condition = "TaskA",
    hrf_estimate = list(NULL),
    taps = list(1:6)
  )
  X_null_hrf <- ndx_build_design_matrix(estimated_hrfs = hrf_null_estimate, events = events_single_dm, motion_params = NULL, rpca_components = NULL,
                                       spectral_sines = NULL, run_idx = run_idx_single_dm, TR = TR_test_dm, poly_degree = -1, verbose = FALSE)
  expect_false(any(startsWith(colnames(X_null_hrf), "task_")))
  expect_equal(ncol(X_null_hrf), 1) # Intercept
  
  # Case 4: HRF estimate is non-numeric
  hrf_non_numeric <- tibble::tibble(
    condition = "TaskA",
    hrf_estimate = list("not_numeric"),
    taps = list(1:6)
  )
  X_non_numeric_hrf <- ndx_build_design_matrix(estimated_hrfs = hrf_non_numeric, events = events_single_dm, motion_params = NULL, rpca_components = NULL,
                                              spectral_sines = NULL, run_idx = run_idx_single_dm, TR = TR_test_dm, poly_degree = -1, verbose = FALSE)
  expect_false(any(startsWith(colnames(X_non_numeric_hrf), "task_")))
  expect_equal(ncol(X_non_numeric_hrf), 1) # Intercept
})

test_that("drop_zero_variance option removes constant regressors", {
  const_rpca <- matrix(1, nrow = total_timepoints_dm, ncol = 1)
  X_drop <- ndx_build_design_matrix(
    estimated_hrfs = estimated_hrfs_dm,
    events = events_dm,
    motion_params = motion_params_dm,
    rpca_components = const_rpca,
    spectral_sines = NULL,
    run_idx = run_idx_dm,
    TR = TR_test_dm,
    poly_degree = 0,
    verbose = FALSE,
    drop_zero_variance = TRUE
  )
  expect_false("rpca_comp_1" %in% colnames(X_drop))
  
  # Test that intercept columns are NOT removed even if constant
  X_keep_intercept <- ndx_build_design_matrix(
    estimated_hrfs = NULL,
    events = events_dm,
    motion_params = NULL,
    rpca_components = const_rpca,
    spectral_sines = NULL,
    run_idx = run_idx_dm,
    TR = TR_test_dm,
    poly_degree = 0,
    verbose = FALSE,
    drop_zero_variance = TRUE
  )
  expect_true("poly0" %in% colnames(X_keep_intercept))
  expect_true("run_intercept_2" %in% colnames(X_keep_intercept))
  expect_false("rpca_comp_1" %in% colnames(X_keep_intercept))
})

test_that("non-sequential run_idx are mapped correctly", {
  run_idx_ns <- rep(c(10, 20), each = n_time_per_run_dm)
  events_ns <- events_dm
  events_ns$blockids <- c(10, 10, 20, 20)

  # Test the internal validation function directly with matching events
  info <- ndx:::.ndx_validate_design_inputs(run_idx_ns, motion_params_dm,
                                            rpca_comps_dm, spectral_sines_dm, events_ns)
  expect_equal(unique(info$run_idx_mapped), c(1, 2))
  expect_equal(info$run_lengths, c(n_time_per_run_dm, n_time_per_run_dm))

  X_ns <- ndx_build_design_matrix(
    estimated_hrfs = estimated_hrfs_dm,
    events = events_ns,
    motion_params = motion_params_dm,
    rpca_components = rpca_comps_dm,
    spectral_sines = spectral_sines_dm,
    run_idx = run_idx_ns,
    TR = TR_test_dm,
    poly_degree = 1,
    verbose = FALSE
  )
  expect_true(is.matrix(X_ns))
  expect_equal(nrow(X_ns), total_timepoints_dm)
})

test_that("ndx_build_design_matrix errors when events blockids are invalid", {
  events_bad <- events_dm
  events_bad$blockids[1] <- 99
  expect_error(
    ndx_build_design_matrix(
      estimated_hrfs = estimated_hrfs_dm,
      events = events_bad,
      motion_params = motion_params_dm,
      rpca_components = rpca_comps_dm,
      spectral_sines = spectral_sines_dm,
      run_idx = run_idx_dm,
      TR = TR_test_dm,
      poly_degree = 1,
      verbose = FALSE
    ),
    "events\\$blockids"
  )
})

test_that("non-sequential run_idx are mapped for events", {
  run_idx_nonseq <- rep(c(10, 20), each = n_time_per_run_dm)
  events_nonseq <- events_dm
  events_nonseq$blockids <- c(10, 10, 20, 20)

  X_nonseq <- ndx_build_design_matrix(
    estimated_hrfs = estimated_hrfs_dm,
    events = events_nonseq,
    motion_params = motion_params_dm,
    rpca_components = rpca_comps_dm,
    spectral_sines = spectral_sines_dm,
    run_idx = run_idx_nonseq,
    TR = TR_test_dm,
    poly_degree = 0,
    verbose = FALSE
  )
  expect_true(is.matrix(X_nonseq))
  expect_equal(nrow(X_nonseq), total_timepoints_dm)
  # Note: The run intercept column name should reflect the mapped run ID (2), not the original ID (20)
  expect_true(any(grepl("run_intercept_2", colnames(X_nonseq))))
})
                         
test_that("near-constant regressors are removed when variance below epsilon", {
  set.seed(123)
  near_const_rpca <- matrix(1 + rnorm(total_timepoints_dm, sd = 1e-10), ncol = 1)
  X_drop_near <- ndx_build_design_matrix(
    estimated_hrfs = estimated_hrfs_dm,
    events = events_dm,
    motion_params = motion_params_dm,
    rpca_components = near_const_rpca,
    spectral_sines = NULL,
    run_idx = run_idx_dm,
    TR = TR_test_dm,
    poly_degree = 0,
    verbose = FALSE,
    drop_zero_variance = TRUE,
    zero_var_epsilon = 1e-8
  )
  expect_false("rpca_comp_1" %in% colnames(X_drop_near))
})

test_that("estimated_hrfs validation works correctly", {
  # Test missing required columns
  bad_hrfs <- tibble::tibble(
    condition = "TaskA",
    # missing hrf_estimate column
    taps = list(1:6)
  )
  expect_error(
    ndx_build_design_matrix(
      estimated_hrfs = bad_hrfs,
      events = events_dm,
      motion_params = NULL,
      rpca_components = NULL,
      spectral_sines = NULL,
      run_idx = run_idx_dm,
      TR = TR_test_dm,
      poly_degree = 0,
      verbose = FALSE
    ),
    "estimated_hrfs tibble must contain 'condition' and 'hrf_estimate' columns"
  )
  
  # Test non-tibble input
  bad_hrfs_df <- data.frame(
    condition = "TaskA",
    hrf_estimate = I(list(rnorm(6)))
  )
  expect_no_error({
    X_df_hrfs <- ndx_build_design_matrix(
      estimated_hrfs = bad_hrfs_df,
      events = events_dm,
      motion_params = NULL,
      rpca_components = NULL,
      spectral_sines = NULL,
      run_idx = run_idx_dm,
      TR = TR_test_dm,
      poly_degree = 0,
      verbose = FALSE
    )
  })
  # Should still work since it's converted to tibble internally
  expect_true(is.matrix(X_df_hrfs))
})

test_that("events validation works correctly", {
  # Test missing blockids column
  bad_events <- data.frame(
    onsets = c(5, 15) * TR_test_dm,
    condition = c("TaskA", "TaskB")
    # missing blockids column
  )
  expect_error(
    ndx_build_design_matrix(
      estimated_hrfs = estimated_hrfs_dm,
      events = bad_events,
      motion_params = NULL,
      rpca_components = NULL,
      spectral_sines = NULL,
      run_idx = run_idx_dm,
      TR = TR_test_dm,
      poly_degree = 0,
      verbose = FALSE
    ),
    "events must be a data.frame containing a 'blockids' column for validation"
  )
  
  # Test non-integer blockids
  bad_events_float <- events_dm
  bad_events_float$blockids[1] <- 1.5
  expect_error(
    ndx_build_design_matrix(
      estimated_hrfs = estimated_hrfs_dm,
      events = bad_events_float,
      motion_params = NULL,
      rpca_components = NULL,
      spectral_sines = NULL,
      run_idx = run_idx_dm,
      TR = TR_test_dm,
      poly_degree = 0,
      verbose = FALSE
    ),
    "events\\$blockids contains values not present in run_idx"
  )
  
  # Test NA blockids
  bad_events_na <- events_dm
  bad_events_na$blockids[1] <- NA
  expect_error(
    ndx_build_design_matrix(
      estimated_hrfs = estimated_hrfs_dm,
      events = bad_events_na,
      motion_params = NULL,
      rpca_components = NULL,
      spectral_sines = NULL,
      run_idx = run_idx_dm,
      TR = TR_test_dm,
      poly_degree = 0,
      verbose = FALSE
    ),
    "events\\$blockids contains values not present in run_idx"
  )
})

test_that("task regressor generation handles onset timing correctly", {
  # Test with events at different time points
  events_timing <- data.frame(
    onsets = c(0, 10, 20, 30) * TR_test_dm, # Events at different times
    durations = rep(2 * TR_test_dm, 4),
    condition = rep("TaskA", 4),
    blockids = c(1, 1, 2, 2)
  )
  
  hrf_single <- create_mock_hrfs(conditions = "TaskA", taps_per_hrf = 4, TR = TR_test_dm)
  
  X_timing <- ndx_build_design_matrix(
    estimated_hrfs = hrf_single,
    events = events_timing,
    motion_params = NULL,
    rpca_components = NULL,
    spectral_sines = NULL,
    run_idx = run_idx_dm,
    TR = TR_test_dm,
    poly_degree = -1,
    verbose = FALSE
  )
  
  expect_true(is.matrix(X_timing))
  expect_equal(ncol(X_timing), 3) # task_TaskA + run_intercept_1 + run_intercept_2
  expect_true("task_TaskA" %in% colnames(X_timing))
  
  # Check that the task regressor has non-zero values at expected time points
  task_col <- X_timing[, "task_TaskA"]
  expect_true(any(task_col != 0)) # Should have some non-zero values
})

test_that("column naming works correctly with special characters", {
  # Test with condition names that need make.names() processing
  special_hrfs <- tibble::tibble(
    condition = c("Task-A", "Task B", "123Task"),
    hrf_estimate = list(rnorm(4), rnorm(4), rnorm(4)),
    taps = list(1:4, 1:4, 1:4)
  )
  
  events_special <- data.frame(
    onsets = c(5, 15, 25) * TR_test_dm,
    durations = rep(2 * TR_test_dm, 3),
    condition = c("Task-A", "Task B", "123Task"),
    blockids = c(1, 1, 2)
  )
  
  X_special <- ndx_build_design_matrix(
    estimated_hrfs = special_hrfs,
    events = events_special,
    motion_params = NULL,
    rpca_components = NULL,
    spectral_sines = NULL,
    run_idx = run_idx_dm,
    TR = TR_test_dm,
    poly_degree = 0,
    verbose = FALSE
  )
  
  expect_true(is.matrix(X_special))
  # Check that column names are valid R names
  expect_true(all(make.names(colnames(X_special)) == colnames(X_special)))
  expect_true(any(grepl("task_Task.A", colnames(X_special))))
  expect_true(any(grepl("task_Task.B", colnames(X_special))))
  expect_true(any(grepl("task_X123Task", colnames(X_special))))
})
