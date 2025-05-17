#' Perform an Initial GLM to Generate Residuals and Baseline Variance
#'
#' This function performs an initial GLM fit to the fMRI data to obtain
#' residuals that will serve as input for subsequent denoising steps (e.g., RPCA).
#' It also calculates the variance of residuals from a simpler task-only model,
#' which is used as the denominator in the Denoising Efficacy Score (DES).
#'
#' @param Y_fmri Matrix of fMRI data (timepoints x voxels). Assumed to be
#'   concatenated across runs if multiple runs are present.
#' @param events A data frame describing experimental events. Must contain columns:
#'   - `onsets`: Numeric, onset time of the event in seconds, relative to the start of its run.
#'   - `durations`: Numeric, duration of the event in seconds.
#'   - `condition`: Character or factor, identifying the event type/condition.
#'   - `blockids`: Numeric, 1-indexed, identifying the run/block each event belongs to.
#'                 This must correspond to the runs in `run_idx` and `Y_fmri`.
#' @param motion_params Matrix or data frame of motion parameters (timepoints x n_regressors).
#'   Must have the same number of rows as `Y_fmri`.
#' @param run_idx Numeric vector indicating run membership for each timepoint in `Y_fmri`.
#'   (e.g., c(rep(1,100), rep(2,100)) for two runs of 100 timepoints each).
#' @param TR Numeric, repetition time in seconds.
#' @return A list containing:
#'   - `Y_residuals_current`: Matrix of residuals after the full initial GLM (task + motion + polynomials).
#'   - `VAR_BASELINE_FOR_DES`: Numeric, variance of residuals from a task-only GLM (task + run intercepts).
#' @examples
#' \dontrun{
#' # Example Usage (requires fmrireg and appropriate data)
#' # Define some hypothetical data
#' n_time_per_run <- 100
#' n_runs <- 2
#' n_voxels <- 10
#' TR <- 2.0
#' total_timepoints <- n_time_per_run * n_runs
#'
#' Y_fmri_example <- matrix(rnorm(total_timepoints * n_voxels), nrow = total_timepoints, ncol = n_voxels)
#' run_idx_example <- rep(1:n_runs, each = n_time_per_run)
#' motion_params_example <- matrix(rnorm(total_timepoints * 6), nrow = total_timepoints, ncol = 6)
#'
#' # Simple events table
#' events_example <- data.frame(
#'   onsets = c(10, 30, 10, 30),
#'   durations = c(5, 5, 5, 5),
#'   condition = rep(c("TaskA", "TaskB"), 2),
#'   blockids = c(1, 1, 2, 2)
#' )
#'
#' initial_glm_output <- ndx_initial_glm(Y_fmri_example, events_example,
#'                                     motion_params_example, run_idx_example, TR)
#' head(initial_glm_output$Y_residuals_current)
#' print(initial_glm_output$VAR_BASELINE_FOR_DES)
#' }
#' @import fmrireg
#' @import stats
#' @export
ndx_initial_glm <- function(Y_fmri, events, motion_params, run_idx, TR) {

  # Validate inputs
  if (!is.matrix(Y_fmri)) {
    stop("Y_fmri must be a matrix (timepoints x voxels).")
  }
  n_timepoints <- nrow(Y_fmri)
  if (n_timepoints == 0) stop("Y_fmri has zero timepoints.")
  if (ncol(Y_fmri) == 0) stop("Y_fmri has zero voxels.")

  if (!is.data.frame(events)) stop("events must be a data frame.")
  required_event_cols <- c("onsets", "durations", "condition", "blockids")
  if (!all(required_event_cols %in% names(events))) {
    stop(paste("events data frame must contain columns:", paste(required_event_cols, collapse=", ")))
  }

  if (!is.matrix(motion_params) && !is.data.frame(motion_params)) {
    stop("motion_params must be a matrix or data frame.")
  }
  if (nrow(motion_params) != n_timepoints) {
    stop("Number of rows in motion_params must match Y_fmri.")
  }

  if (!is.numeric(run_idx) || length(run_idx) != n_timepoints) {
    stop("run_idx must be a numeric vector with length matching nrow(Y_fmri).")
  }
  if (!is.numeric(TR) || length(TR) != 1 || TR <= 0) {
    stop("TR must be a single positive number.")
  }

  # 1. Create Sampling Frame
  # table() on run_idx gives counts per run, which are the blocklens
  # Ensure runs are ordered if run_idx is not already e.g. 1,1,1,2,2,2
  unique_runs <- unique(run_idx)
  if (any(diff(unique_runs) < 0) && length(unique_runs) >1 ) { # check if not sorted if multiple runs
      # This ensures blocklens are in the order of appearance of runs
      run_lengths <- as.numeric(table(factor(run_idx, levels=unique(run_idx))))
  } else {
      run_lengths <- as.numeric(table(run_idx))
  }
  
  if(any(run_lengths == 0)) stop("Some runs specified in run_idx have zero length.")

  sf <- fmrireg::sampling_frame(blocklens = run_lengths, TR = TR)

  # 2. Task-Only Model (for VAR_BASELINE_FOR_DES)
  task_conditions <- unique(as.character(events$condition))
  
  # Create a copy of events to modify for the model
  events_for_model_wide <- events
  if (!is.factor(events_for_model_wide$condition)) {
    events_for_model_wide$condition <- factor(events_for_model_wide$condition)
  }

  if (length(task_conditions) == 0) {
    warning("No task conditions found in 'events$condition'. VAR_BASELINE_FOR_DES will be based on a baseline-only model.")
    task_formula_str <- "~ 0" 
  } else {
    formula_terms <- character(length(task_conditions))
    for (i in seq_along(task_conditions)) {
      cond_name <- task_conditions[i]
      # Create a syntactically valid column name for the dummy variable
      dummy_col_name <- make.names(cond_name) 
      
      # Add a dummy coded column for the current condition to the events data frame
      # This column will be TRUE/1 for rows matching the current condition, FALSE/0 otherwise.
      # This is NOT what fmrireg::hrf expects. hrf(var) expects 'var' to define onsets for that event type.
      # The original 'events' table is long format. fmrireg should use onsets/durations from it, filtered by condition.

      # Let's revert to the idea that hrf(condition_factor_column_name, ...) should work if fmrireg is robust.
      # The previous errors might have been due to something else if this structure is standard for fmrireg.
      # The most standard fmrireg usage for multiple conditions from a single factor column in `data`
      # when onsets/durations are also columns in `data` is indeed `~ hrf(condition_factor_column_name, ...)`
      # where `condition_factor_column_name` is the actual name of the column in the `data` data.frame.
    }
    # This simplified approach should be the one fmrireg handles if it supports long-format event tables well.
    task_formula_str <- "onsets ~ fmrireg::hrf(condition, basis=\"spmg1\")"
  }
  
  task_formula <- as.formula(task_formula_str)

  # Pass the events table (with condition as factor).
  # fmrireg::event_model will use the 'onsets' column from LHS of formula,
  # and 'condition' for grouping from RHS, all from the 'data' argument.
  # The 'block' argument aligns these to runs in the sampling_frame.
  em_task_only <- fmrireg::event_model(task_formula, 
                                     data = events_for_model_wide, 
                                     block = events_for_model_wide$blockids, 
                                     sampling_frame = sf)
  bm_task_only <- fmrireg::baseline_model(basis = "constant", intercept = "runwise", sframe = sf)
  model_task_only <- fmrireg::fmri_model(event_model = em_task_only, baseline_model = bm_task_only)
  
  X_task_only <- tryCatch({
    fmrireg::design_matrix(model_task_only)
  }, error = function(e) {
    stop(paste("Error creating task-only design matrix:", e$message))
  })
  
  if (ncol(X_task_only) == 0 && length(task_conditions) > 0) {
      warning("Task-only design matrix has zero columns despite conditions being present. Check event timings and model spec.")
  } else if (ncol(X_task_only) == 0 && length(task_conditions) == 0) {
      warning("Task-only design matrix (baseline only) has zero columns. This is unexpected.")
  }

  resid_task_only <- calculate_residuals_ols(Y_fmri, as.matrix(X_task_only))
  VAR_BASELINE_FOR_DES <- stats::var(as.vector(resid_task_only))

  # 3. Pass 0 Model (for Y_residuals_current)
  # Event model (em_task_only) is the same.
  # Baseline model includes motion and polynomials.
  # Ensure motion_params is a matrix for nuisance_list
  mat_motion_params <- as.matrix(motion_params)
  # Column names for motion parameters are tricky if they are to be split into a list for nuisance_list.
  # fmrireg might expect consistent naming or just a list of matrices.
  # For now, we'll name the concatenated matrix, then split.
  if (ncol(mat_motion_params) > 0) {
      colnames(mat_motion_params) <- paste0("motion_col_", 1:ncol(mat_motion_params))
  }

  # Prepare nuisance_list based on number of runs
  n_runs_from_sf <- length(sf$blocklens) # Get number of runs from sampling_frame
  nuisance_arg_for_bm_pass0 <- NULL
  if (ncol(mat_motion_params) > 0) { # Only add nuisance if there are motion params
    if (n_runs_from_sf > 1) {
      # Split mat_motion_params by run_idx into a list of matrices
      # Each element of the list will be the motion regressors for that run
      split_motion_params <- split.data.frame(mat_motion_params, run_idx)
      # Convert data frames in list to matrices
      nuisance_arg_for_bm_pass0 <- lapply(split_motion_params, as.matrix)
      # fmrireg might expect this to be a named list if multiple nuisance types, 
      # or an unnamed list if it's just one type of nuisance split by block.
      # The error message suggests the *number* of elements in nuisance_list must match blocks.
      # If nuisance_list itself should contain one element per block for a *given* nuisance type,
      # then this structure (a list of matrices) is what fmrireg expects for `list(motion_effects=LIST_OF_MATRICES_PER_RUN)`.
      # Or, does it expect nuisance_list = list(run1_motion=mat1, run2_motion=mat2)?
      # Let's try making nuisance_list a list where each element is a nuisance matrix for a run.
      # The error says "number of `nuisance_list` elements must match number of blocks"
      # This implies nuisance_list = list(run1_nuis, run2_nuis, ...)
      # where each run_i_nuis is a matrix. Let's test this interpretation.
      # The elements of this list should probably be matrices of nuisance for each run.
      # So, if `mat_motion_params` is the only nuisance, we pass a list of these split matrices.
    } else { 
      # Single run, nuisance_list has one element which is the full (single-run) motion matrix
      nuisance_arg_for_bm_pass0 <- list(motion_effects_run1 = mat_motion_params)
    }
  }

  bm_pass0 <- fmrireg::baseline_model(basis = "poly", degree = 1L, intercept = "runwise",
                                      sframe = sf, nuisance_list = nuisance_arg_for_bm_pass0)

  model_pass0 <- fmrireg::fmri_model(event_model = em_task_only, baseline_model = bm_pass0)
  
  X_pass0 <- tryCatch({
    fmrireg::design_matrix(model_pass0)
  }, error = function(e) {
    stop(paste("Error creating Pass-0 design matrix:", e$message))
  })

  if (ncol(X_pass0) == 0) {
      stop("Pass-0 design matrix has zero columns. This indicates a serious issue with model specification.")
  }
  
  Y_residuals_current <- calculate_residuals_ols(Y_fmri, as.matrix(X_pass0))

  return(list(
    Y_residuals_current = Y_residuals_current,
    VAR_BASELINE_FOR_DES = VAR_BASELINE_FOR_DES
  ))
} 