#' Estimate Initial FIR HRFs using Fused Lasso
#'
#' Estimates a single "global" or "good-voxel-average" FIR HRF per condition
#' using `glmgen::fusedlasso`. The HRF is estimated from the robust mean time course
#' of "good" voxels, identified by R-squared from Pass 0 residuals.
#' Lambda parameters for fusedlasso are chosen via block-wise k-fold cross-validation.
#'
#' @param Y_fmri Matrix of fMRI data (timepoints x voxels), concatenated across runs.
#' @param pass0_residuals Matrix of residuals from `ndx_initial_glm` (timepoints x voxels).
#' @param events A data frame describing experimental events. Must contain columns:
#'   `onsets`, `durations`, `condition`, `blockids`.
#' @param run_idx Numeric vector indicating run membership for each timepoint in `Y_fmri`.
#' @param TR Numeric, repetition time in seconds.
#' @param spike_TR_mask Optional. A logical vector of length `nrow(Y_fmri)` where TRUE
#'   indicates a TR to be excluded from HRF estimation (e.g., due to spikes).
#'   If NULL, all TRs are considered valid.
#' @param user_options A list of user-configurable options:
#'   - `hrf_fir_taps` (integer): Number of FIR basis functions (e.g., 20).
#'   - `hrf_fir_span_seconds` (numeric): Total duration the FIR model covers (e.g., 24s).
#'   - `good_voxel_R2_threshold` (numeric): R-squared threshold for selecting good voxels (e.g., 0.02-0.06).
#'   - `cv_folds` (integer): Number of folds for K-fold cross-validation (e.g., 5).
#'   - `lambda1_grid` (numeric vector): Grid of values for fused penalty `lambda` (L1 on differences).
#'   - `lambda2_grid` (numeric vector): Grid of values for ridge penalty `gamma` (L2 on coefficients).
#'   - `hrf_min_good_voxels` (integer): Minimum number of good voxels required to proceed (e.g., 100).
#'   - `return_full_model` (logical): If TRUE, include the full glmgen model object in the output. Defaults to FALSE.
#' @return A tibble with columns: `condition` (character), `hrf_estimate` (list of numeric vectors),
#'   `taps` (list of integer vectors). If `user_options$return_full_model` is TRUE, an additional
#'   list-column `glmgen_fit` is included. Returns NULL if critical errors occur.
#' @examples
#' \dontrun{
#' # --- Setup from ndx_initial_glm example --- 
#' n_time_per_run <- 100; n_runs <- 2; n_voxels <- 200; TR <- 2.0
#' total_timepoints <- n_time_per_run * n_runs
#' Y_fmri_example <- matrix(rnorm(total_timepoints * n_voxels), 
#'                          nrow = total_timepoints, ncol = n_voxels)
#' run_idx_example <- rep(1:n_runs, each = n_time_per_run)
#' motion_params_example <- matrix(rnorm(total_timepoints * 6), 
#'                                 nrow = total_timepoints, ncol = 6)
#' events_example <- data.frame(
#'   onsets = c(10, 30, 50, 10, 30, 50),
#'   durations = c(5, 5, 5, 5, 5, 5),
#'   condition = rep(c("TaskA", "TaskB"), each = 3),
#'   blockids = c(1, 1, 1, 2, 2, 2)
#' )
#' # Generate some real signal for TaskA in first 50 voxels
#' sf_hrf <- fmrireg::sampling_frame(blocklens = as.numeric(table(run_idx_example)), TR=TR)
#' design_taskA <- fmrireg::event_model(~ fmrireg::hrf(TaskA, basis="spmg1"), 
#'                                      data=events_example[events_example$condition == "TaskA",], 
#'                                      block=events_example$blockids[events_example$condition == "TaskA"],
#'                                      sampling_frame=sf_hrf)
#' X_taskA_signal <- fmrireg::design_matrix(design_taskA)
#' if (ncol(X_taskA_signal) > 0) { 
#'    Y_fmri_example[,1:50] <- Y_fmri_example[,1:50] + X_taskA_signal %*% matrix(rnorm(ncol(X_taskA_signal)*50, mean=3, sd=1), ncol=50)
#' }
#' 
#' pass0_out <- ndx_initial_glm(Y_fmri_example, events_example, 
#'                                  motion_params_example, run_idx_example, TR)
#' 
#' # --- Now for ndx_estimate_initial_hrfs --- 
#' user_opts_hrf <- list(
#'   hrf_fir_taps = 12,            
#'   hrf_fir_span_seconds = 24,  
#'   good_voxel_R2_threshold = 0.01,
#'   cv_folds = 2, 
#'   lambda1_grid = 10^seq(-2, 1, length.out = 3),
#'   lambda2_grid = 10^seq(-3, 0, length.out = 3),
#'   hrf_min_good_voxels = 10,
#'   return_full_model = FALSE
#' )
#' 
#' estimated_hrfs_tbl <- ndx_estimate_initial_hrfs(
#'   Y_fmri = Y_fmri_example,
#'   pass0_residuals = pass0_out$Y_residuals_current,
#'   events = events_example,
#'   run_idx = run_idx_example,
#'   TR = TR,
#'   spike_TR_mask = NULL,
#'   user_options = user_opts_hrf
#' )
#' 
#' if (!is.null(estimated_hrfs_tbl)) {
#'   print(estimated_hrfs_tbl)
#'   # Example of unnesting for plotting (requires tidyr and ggplot2)
#'   # if (requireNamespace("tidyr", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)) {
#'   #   plot_df <- tidyr::unnest(estimated_hrfs_tbl, cols = c(hrf_estimate, taps))
#'   #   print(ggplot2::ggplot(plot_df, ggplot2::aes(x = taps * TR, y = hrf_estimate)) +
#'   #         ggplot2::geom_line() + ggplot2::geom_point() +
#'   #         ggplot2::facet_wrap(~condition) + ggplot2::labs(x="Time (s)", y="Amplitude"))
#'   # }
#' }
#' }
#' @import fmrireg
#' @import glmgen
#' @importFrom genlasso fusedlasso
#' @import stats
#' @import matrixStats
#' @import tibble
#' @export
ndx_estimate_initial_hrfs <- function(Y_fmri, pass0_residuals, events,
                                        run_idx, TR,
                                        spike_TR_mask = NULL,
                                        user_options) {

  validated_inputs <- validate_hrf_inputs(Y_fmri, pass0_residuals, events, run_idx, TR, 
                                            spike_TR_mask, user_options)
  user_options <- validated_inputs$user_options
  # n_timepoints <- validated_inputs$n_timepoints 
  # spike_TR_mask in this scope is the validated one from validated_inputs

  response_data <- prepare_hrf_response_data(Y_fmri, pass0_residuals, run_idx, 
                                             validated_inputs$validated_spike_TR_mask, 
                                             user_options)
  
  if (is.null(response_data)) {
    return(NULL) 
  }
  
  ybar_clean <- response_data$ybar_clean
  block_ids_for_cv <- response_data$block_ids_for_cv
  valid_TRs_for_hrf_estimation_mask <- response_data$valid_TRs_for_hrf_estimation_mask
  
  results_list_collector <- list() # Changed name to avoid clash with internal var name
  conditions <- unique(as.character(events$condition))
  
  unique_runs_overall <- unique(run_idx)
  run_lengths_overall <- as.numeric(table(factor(run_idx, levels = unique_runs_overall)))
  sf_overall <- fmrireg::sampling_frame(blocklens = run_lengths_overall, TR = TR)
  
  for (cond_name_iter in conditions) { # Renamed loop variable
    message(paste("Processing HRF for condition:", cond_name_iter))
    events_for_this_condition <- events[events$condition == cond_name_iter, , drop = FALSE]
    
    # Call the new per-condition estimation function
    condition_hrf_result <- estimate_hrf_for_condition(
      condition_name = cond_name_iter,
      events_for_condition = events_for_this_condition,
      ybar_clean = ybar_clean,
      block_ids_for_cv = block_ids_for_cv,
      overall_sampling_frame = sf_overall,
      valid_TRs_mask = valid_TRs_for_hrf_estimation_mask, # pass the mask directly
      TR = TR,
      user_options = user_options
    )
    results_list_collector[[cond_name_iter]] <- condition_hrf_result
  }
  
  if (length(results_list_collector) == 0) return(NULL) # Should be caught if conditions is empty
  
  # Convert list of lists to a tibble
  # Filter out any top-level NULLs if a condition somehow failed to even produce a list structure (unlikely with current setup)
  valid_results <- Filter(Negate(is.null), results_list_collector)
  if(length(valid_results) == 0 && length(conditions) > 0) {
      warning("No valid HRF results were generated for any condition.")
      return(NULL)
  }
  if(length(valid_results) == 0 && length(conditions) == 0) {
      message("No conditions found in event data to estimate HRFs for.")
      return(NULL)
  }

  # Reconstruct based on the structure returned by estimate_hrf_for_condition
  conditions_out <- sapply(valid_results, function(x) x$condition)
  hrf_estimates_out <- lapply(valid_results, function(x) x$hrf_estimate)
  taps_out <- lapply(valid_results, function(x) x$taps)
  
  output_tbl <- tibble::tibble(
      condition = conditions_out,
      hrf_estimate = hrf_estimates_out,
      taps = taps_out
  )
  
  if (user_options$return_full_model) {
      glmgen_fits_out <- lapply(valid_results, function(x) x$glmgen_fit)
      # Ensure this column is added only if there are valid results
      if (nrow(output_tbl) > 0) {
        output_tbl$glmgen_fit <- glmgen_fits_out
      } else if (length(glmgen_fits_out) > 0) {
        # This case (no rows in output_tbl but fits exist) should not happen if logic is correct
        # but as a fallback, create a tibble for fits if conditions_out was empty but fits were somehow produced
        output_tbl <- tibble::tibble(condition=NA_character_, hrf_estimate=list(NULL), taps=list(NULL), glmgen_fit = glmgen_fits_out)[-1,]
        warning("HRF fitting produced glmgen_fits but no primary HRF data; check for inconsistencies.")
      }
  }
  
  if (nrow(output_tbl) == 0 && length(conditions) > 0) {
      warning("HRF estimation resulted in an empty table despite conditions being present. Check warnings for individual conditions.")
      return(NULL)
  }

  return(output_tbl)
}

#' Generate FIR Design Matrix for a Specific Condition
#' @keywords internal
get_fir_design_matrix_for_condition <- function(condition_name, events_df,
                                                sampling_frame, fir_taps, TR) {
  if (nrow(events_df) == 0) {
    return(matrix(0, nrow = sampling_frame$total_samples, ncol = fir_taps))
  }
  
  safe_cond_name <- make.names(condition_name, unique = FALSE)
  fir_span <- fir_taps * TR 
  
  if (safe_cond_name != condition_name) {
    formula_term <- sprintf("fmrireg::hrf(`%s`, basis=\"tent\", nbasis=%d, span=%.2f)", 
                            condition_name, fir_taps, fir_span)
  } else {
    formula_term <- sprintf("fmrireg::hrf(%s, basis=\"tent\", nbasis=%d, span=%.2f)", 
                            condition_name, fir_taps, fir_span)
  }
  
  model_formula_str <- paste("~", formula_term)
  model_formula <- as.formula(model_formula_str)
  
  ev_model <- fmrireg::event_model(model_formula,
                                 data = events_df, 
                                 block = events_df$blockids,
                                 sampling_frame = sampling_frame,
                                 drop_empty = FALSE) 
  
  X_fir <- fmrireg::design_matrix(ev_model)
  
  if (ncol(X_fir) != fir_taps) {
      if (ncol(X_fir) == fir_taps + 1 && all(X_fir[,1] == 1)) { 
          warning(sprintf("FIR design for %s included an intercept; removing it.", condition_name))
          X_fir <- X_fir[, -1, drop=FALSE]
      } else if (ncol(X_fir) == 0 && fir_taps > 0 && nrow(events_df) > 0) {
          warning(sprintf("FIR design for %s (with %d events) resulted in 0 columns. Expected %d. Check event timings relative to scan length and FIR span. Returning zeros.", 
                          condition_name, nrow(events_df), fir_taps))
          return(matrix(0, nrow = sampling_frame$total_samples, ncol = fir_taps))
      } else if (ncol(X_fir) != fir_taps && nrow(events_df) > 0) {
        stop(sprintf("FIR design for %s has %d columns, expected %d. FIR span %.2f. Events: %d. Check event data relative to FIR parameters and scan times.", 
                     condition_name, ncol(X_fir), fir_taps, fir_span, nrow(events_df)))
      } else if (ncol(X_fir) == 0 && fir_taps > 0 && nrow(events_df) == 0) {
         return(matrix(0, nrow = sampling_frame$total_samples, ncol = fir_taps))
      }
  }
  
  return(X_fir)
}

#' Basic Cone Projection Heuristic for FIR Coefficients
#' @keywords internal
project_cone_heuristic <- function(fir_coeffs) {
  if (length(fir_coeffs) == 0) return(numeric(0))
  
  coeffs_adj <- fir_coeffs
  if (length(coeffs_adj) > 0) { 
    coeffs_adj <- coeffs_adj - coeffs_adj[1]
  }
  
  if (length(coeffs_adj) >= 2 && coeffs_adj[2] < 0) {
    coeffs_adj[2] <- 0
  }
  
  # If HRF became flat (all zeros) after initial adjustment, give it a small bump.
  if (sum(abs(coeffs_adj)) < 1e-9 && length(coeffs_adj) > 0) {
      # Place a small positive value at a plausible early peak (e.g., 1/3rd of the way through taps)
      peak_tap_idx <- max(1, floor(length(coeffs_adj) / 3))
      coeffs_adj[peak_tap_idx] <- 1e-6
  }
  
  return(coeffs_adj)
}

#' Block-wise K-fold Cross-Validation for Fused Lasso (glmgen)
#' @keywords internal
cv_fusedlasso <- function(y, X, lambda_grid, gamma_grid, k_folds, block_ids) {
  if (nrow(X) != length(y)) stop("nrow(X) must equal length(y) for CV.")
  if (length(block_ids) != length(y)) stop("Length of block_ids must match length(y).")
  
  p_coeffs <- ncol(X)
  D_1d <- NULL
  if (p_coeffs > 1) {
    D_1d <- diff(diag(p_coeffs), differences = 1)
  } else if (p_coeffs == 1) {
    D_1d <- diff(diag(1), differences = 1) # 0-row matrix, results in no fusion penalty
  }

  unique_blocks <- unique(block_ids)
  n_unique_blocks <- length(unique_blocks)
  
  if (n_unique_blocks == 0) stop("No blocks provided for CV.")
  if (n_unique_blocks < 2 && k_folds > 1) {
      warning("Only 1 unique block for CV. Cannot perform k-fold CV with k > 1. Setting k_folds = 1 (no CV, using first lambda/gamma).")
      k_folds <- 1 
  }
  if (k_folds > n_unique_blocks) {
    warning(paste("k_folds (", k_folds, ") is greater than the number of unique blocks (", n_unique_blocks,
                  "). Setting k_folds to ", n_unique_blocks, " (leave-one-block-out CV).", sep=""))
    k_folds <- n_unique_blocks
  }
  if (nrow(X) < k_folds && k_folds > 1) { 
      warning(paste("Number of observations (", nrow(X), ") is less than k_folds (", k_folds, ") after block considerations. Adjusting k_folds."))
      k_folds <- max(1, min(k_folds, nrow(X))) 
  }

  fold_assignment_for_blocks <- sample(rep(1:k_folds, length.out = n_unique_blocks))
  observation_fold_indices <- fold_assignment_for_blocks[match(block_ids, unique_blocks)]
  
  cv_errors <- array(NA, dim = c(length(lambda_grid), length(gamma_grid)))
  dimnames(cv_errors) <- list(lambda = signif(lambda_grid,3), gamma = signif(gamma_grid,3))

  for (i_g in seq_along(gamma_grid)) { # Outer loop for gamma (L2 penalty)
    gamma_val <- gamma_grid[i_g]
    fold_mse_for_gamma <- matrix(NA, nrow = k_folds, ncol = length(lambda_grid))

    # Pre-calculate the fusedlasso path for this gamma_val if k_folds > 1
    # For k_folds = 1, we do it inside the lambda loop as it's just one fit.

    for (k in 1:k_folds) {
      test_idx <- which(observation_fold_indices == k)
      train_idx <- which(observation_fold_indices != k)

      if (k_folds == 1) { # Special case: no actual CV, fit on full data for each lambda/gamma
          # This logic might need refinement if k_folds=1 means predict on same data
          # For now, assume it means fit on y, X and then evaluate for each lambda
          y_train <- y; X_train <- X
          y_test <- y; X_test <- X # Evaluate on same data if k_folds = 1
      } else {
          if (length(train_idx) == 0 || length(test_idx) == 0) {
              fold_mse_for_gamma[k, ] <- Inf 
              next
          }
          y_train <- y[train_idx]
          X_train <- X[train_idx, , drop = FALSE]
          y_test <- y[test_idx]
          X_test <- X[test_idx, , drop = FALSE]
      }
      
      if (nrow(X_train) < ncol(X_train) || nrow(X_train) < 2) { 
          fold_mse_for_gamma[k, ] <- Inf
          next
      }

      fit_path <- tryCatch({
        genlasso::fusedlasso(y = y_train, X = X_train, D = D_1d, gamma = gamma_val)
      }, error = function(e) {
        warning(sprintf("genlasso::fusedlasso failed for gamma=%.4f on fold %d: %s. Skipping.", gamma_val, k, e$message))
        return(NULL)
      })

      if (is.null(fit_path) || !inherits(fit_path, "genlasso")) {
        fold_mse_for_gamma[k, ] <- Inf
        next
      }

      for (i_l in seq_along(lambda_grid)) {
        lambda_val <- lambda_grid[i_l]
        beta_coeffs <- tryCatch({
          stats::coef(fit_path, lambda = lambda_val)$beta
        }, error = function(e) { NULL })

        if (is.null(beta_coeffs) || length(beta_coeffs) != ncol(X_test)) {
          fold_mse_for_gamma[k, i_l] <- Inf
        } else {
          preds <- X_test %*% beta_coeffs
          fold_mse_for_gamma[k, i_l] <- mean((y_test - preds)^2, na.rm = TRUE)
        }
      }
    } # end k_folds loop
    
    # Average MSE across folds for each lambda, for the current gamma
    for (i_l in seq_along(lambda_grid)) {
        current_lambda_mses <- fold_mse_for_gamma[, i_l]
        mean_mse <- mean(current_lambda_mses[is.finite(current_lambda_mses)], na.rm = TRUE)
        if (!is.finite(mean_mse)) mean_mse <- Inf
        cv_errors[i_l, i_g] <- mean_mse
    }
  } # end gamma_grid loop
  
  if (all(!is.finite(cv_errors))) {
      warning("All CV folds resulted in errors or non-finite MSE for all lambda/gamma combinations. Cannot select optimal parameters.")
      best_l <- lambda_grid[1] 
      best_g <- gamma_grid[1]  
  } else {
      min_error_idx <- arrayInd(which.min(cv_errors), dim(cv_errors))
      best_l <- lambda_grid[min_error_idx[1]]
      best_g <- gamma_grid[min_error_idx[2]]
  }
  
  return(list(best_lambda = best_l, 
              best_gamma = best_g, 
              cv_error_matrix = cv_errors))
}

#' Validate Inputs for HRF Estimation
#'
#' Helper function to validate inputs for `ndx_estimate_initial_hrfs`.
#' Also merges user_options with defaults.
#'
#' @param Y_fmri Matrix of fMRI data.
#' @param pass0_residuals Matrix of residuals.
#' @param events Data frame of experimental events.
#' @param run_idx Numeric vector of run membership.
#' @param TR Numeric, repetition time.
#' @param spike_TR_mask Optional logical vector for spike TRs.
#' @param user_options A list of user-configurable options.
#' @return A list containing validated and potentially augmented `user_options`,
#'   and `n_timepoints` and `validated_spike_TR_mask`.
#' @keywords internal
validate_hrf_inputs <- function(Y_fmri, pass0_residuals, events, run_idx, TR, spike_TR_mask, user_options) {
  if (missing(user_options)) stop("user_options are required for HRF estimation.")
  
  default_opts <- list(
    hrf_fir_taps = 12L,
    hrf_fir_span_seconds = 24,
    good_voxel_R2_threshold = 0.05,
    cv_folds = 5L,
    lambda1_grid = 10^seq(-2, 1, length.out = 5),
    lambda2_grid = 10^seq(-3, 0, length.out = 5),
    hrf_min_good_voxels = 50L,
    return_full_model = FALSE
  )
  user_options <- utils::modifyList(default_opts, user_options)

  if (!is.matrix(Y_fmri) || !is.numeric(Y_fmri)) stop("Y_fmri must be a numeric matrix.")
  n_timepoints <- nrow(Y_fmri)
  if (n_timepoints == 0) stop("Y_fmri has zero timepoints.")
  if (ncol(Y_fmri) == 0) stop("Y_fmri has zero voxels.")

  if (!is.matrix(pass0_residuals) || !is.numeric(pass0_residuals) || 
      nrow(pass0_residuals) != n_timepoints || ncol(pass0_residuals) != ncol(Y_fmri)) {
    stop("pass0_residuals must be a numeric matrix with dimensions matching Y_fmri.")
  }

  if (!is.data.frame(events)) stop("events must be a data frame.")
  required_event_cols <- c("onsets", "durations", "condition", "blockids")
  if (!all(required_event_cols %in% names(events))) {
    stop(paste("events data frame must contain columns:", paste(required_event_cols, collapse=", ")))
  }
  if(!is.numeric(events$onsets) || !is.numeric(events$durations) || !is.numeric(events$blockids)) {
      stop("Event columns 'onsets', 'durations', and 'blockids' must be numeric.")
  }

  if (!is.numeric(run_idx) || length(run_idx) != n_timepoints) {
    stop("run_idx must be a numeric vector with length matching nrow(Y_fmri).")
  }
  if (!is.numeric(TR) || length(TR) != 1 || TR <= 0) {
    stop("TR must be a single positive number.")
  }

  validated_spike_TR_mask <- spike_TR_mask
  if (is.null(validated_spike_TR_mask)) {
    validated_spike_TR_mask <- rep(FALSE, n_timepoints)
  }
  if (length(validated_spike_TR_mask) != n_timepoints || !is.logical(validated_spike_TR_mask)) {
      stop("spike_TR_mask must be a logical vector with length matching n_timepoints in Y_fmri.")
  }
  
  # Validate specific user_options types and values
  if(!is.numeric(user_options$hrf_fir_taps) || user_options$hrf_fir_taps <=0) stop("hrf_fir_taps must be a positive integer.")
  if(!is.numeric(user_options$hrf_fir_span_seconds) || user_options$hrf_fir_span_seconds <=0) stop("hrf_fir_span_seconds must be a positive number.")
  if(!is.numeric(user_options$good_voxel_R2_threshold) || user_options$good_voxel_R2_threshold < 0 || user_options$good_voxel_R2_threshold > 1) stop("good_voxel_R2_threshold must be between 0 and 1.")
  if(!is.numeric(user_options$cv_folds) || user_options$cv_folds <=0) stop("cv_folds must be a positive integer.")
  if(!is.numeric(user_options$lambda1_grid) || length(user_options$lambda1_grid)==0) stop("lambda1_grid must be a non-empty numeric vector.")
  if(!is.numeric(user_options$lambda2_grid) || length(user_options$lambda2_grid)==0) stop("lambda2_grid must be a non-empty numeric vector.")
  if(!is.numeric(user_options$hrf_min_good_voxels) || user_options$hrf_min_good_voxels <0) stop("hrf_min_good_voxels must be a non-negative integer.")
  if(!is.logical(user_options$return_full_model) || length(user_options$return_full_model)!=1) stop("return_full_model must be a single logical value.")

  return(list(user_options = user_options, 
              n_timepoints = n_timepoints, 
              validated_spike_TR_mask = validated_spike_TR_mask))
}

#' Prepare Response Data for HRF Estimation
#'
#' Calculates R2, selects good voxels, and prepares the robust mean time course (`ybar_clean`)
#' from these good voxels, along with associated masks and block IDs for cross-validation.
#'
#' @param Y_fmri Validated fMRI data matrix.
#' @param pass0_residuals Validated Pass 0 residuals matrix.
#' @param run_idx Validated run index vector.
#' @param validated_spike_TR_mask Validated logical mask for spike TRs.
#' @param user_options Validated list of user options, must include 
#'   `good_voxel_R2_threshold` and `hrf_min_good_voxels`.
#' @return A list containing:
#'   - `ybar_clean`: Numeric vector, the robust mean of good voxels on valid TRs.
#'   - `block_ids_for_cv`: Numeric vector, run/block IDs for `ybar_clean`.
#'   - `valid_TRs_for_hrf_estimation_mask`: Logical mask for all original TRs indicating validity for HRF estimation.
#'   - `n_good_voxels`: Integer, count of selected good voxels.
#'   Returns NULL if insufficient good voxels or no valid TRs for estimation.
#' @keywords internal
prepare_hrf_response_data <- function(Y_fmri, pass0_residuals, run_idx, validated_spike_TR_mask, user_options) {
  R2_pass0_vox <- calculate_R2_voxelwise(Y_fmri, pass0_residuals)
  good_voxels_idx <- which(R2_pass0_vox >= user_options$good_voxel_R2_threshold)
  n_good_voxels <- length(good_voxels_idx)

  if (n_good_voxels < user_options$hrf_min_good_voxels) {
    warning(sprintf("Insufficient good voxels found (%d) based on R2 threshold (%.2f). Minimum required: %d. Skipping HRF estimation.",
                    n_good_voxels, user_options$good_voxel_R2_threshold, user_options$hrf_min_good_voxels))
    return(NULL)
  }
  
  Y_good_voxels <- Y_fmri[, good_voxels_idx, drop = FALSE]
  if (ncol(Y_good_voxels) == 0) { 
      warning("No good voxels selected (ncol=0) despite passing minimum count. Check R2 threshold or data. Skipping HRF estimation.")
      return(NULL)
  }

  valid_TRs_for_hrf_estimation_mask <- !validated_spike_TR_mask
  if (sum(valid_TRs_for_hrf_estimation_mask) == 0) {
      warning("No valid TRs available for HRF estimation after applying spike_TR_mask. Skipping HRF estimation.")
      return(NULL)
  }
  
  # Calculate rowMedians on the subset of good voxels, for all original timepoints
  ybar_all_trs <- matrixStats::rowMedians(Y_good_voxels, na.rm = TRUE) 
  # Then select only the valid (non-spike) TRs for ybar_clean
  ybar_clean <- ybar_all_trs[valid_TRs_for_hrf_estimation_mask]
  # Corresponding block IDs for these valid TRs
  block_ids_for_cv <- run_idx[valid_TRs_for_hrf_estimation_mask]
  
  return(list(
    ybar_clean = ybar_clean,
    block_ids_for_cv = block_ids_for_cv,
    valid_TRs_for_hrf_estimation_mask = valid_TRs_for_hrf_estimation_mask,
    n_good_voxels = n_good_voxels
  ))
}

#' Estimate FIR HRF for a Single Condition
#'
#' Performs FIR design matrix generation, cross-validation for fused lasso parameters,
#' final model fitting, and post-processing for a single experimental condition.
#'
#' @param condition_name Character, the name of the condition.
#' @param events_for_condition Data frame of events, already filtered for this specific condition.
#' @param ybar_clean Numeric vector, the robust mean time course of good voxels on valid TRs.
#' @param block_ids_for_cv Numeric vector, run/block IDs for `ybar_clean`.
#' @param overall_sampling_frame An `fmri_sampling_frame` object for all original TRs.
#' @param valid_TRs_mask Logical mask for all original TRs indicating validity for HRF estimation.
#' @param TR Numeric, repetition time.
#' @param user_options List of user options, including `hrf_fir_taps`, `lambda1_grid`, 
#'   `lambda2_grid`, `cv_folds`, and `return_full_model`.
#' @return A list containing:
#'   - `condition`: Character, the condition name.
#'   - `hrf_estimate`: Numeric vector of estimated FIR coefficients, or NULL if estimation fails.
#'   - `taps`: Integer vector of tap indices, or NULL.
#'   - `glmgen_fit`: The `glmgen::fusedlasso` model object if `user_options$return_full_model` is TRUE, else NULL.
#' @keywords internal
estimate_hrf_for_condition <- function(condition_name, events_for_condition,
                                       ybar_clean, block_ids_for_cv,
                                       overall_sampling_frame, valid_TRs_mask, TR,
                                       user_options) {
  
  current_result <- list(condition = condition_name, hrf_estimate = NULL, taps = NULL)
  if (user_options$return_full_model) current_result$glmgen_fit <- NULL

  if (nrow(events_for_condition) == 0) {
      warning(paste("No events found for condition (passed to estimate_hrf_for_condition):", condition_name)) # Should be caught earlier
      return(current_result)
  }
  
  X_fir_cond_all_trs <- tryCatch({
       get_fir_design_matrix_for_condition(
          condition_name = condition_name,
          events_df = events_for_condition, 
          sampling_frame = overall_sampling_frame, 
          fir_taps = user_options$hrf_fir_taps,
          TR = TR 
      )
  }, error = function(e) {
      warning(paste("Error generating FIR design for condition", condition_name, "in estimate_hrf_for_condition:", e$message))
      return(NULL)
  })

  if (is.null(X_fir_cond_all_trs) || ncol(X_fir_cond_all_trs) == 0) {
    warning(paste("FIR design matrix for condition", condition_name, "is empty or NULL in estimate_hrf_for_condition."))
    return(current_result)
  }
  
  # Subset the design matrix to include only non-spike TRs (valid for estimation)
  X_fir_cond_clean <- X_fir_cond_all_trs[valid_TRs_mask, , drop = FALSE]
  
  if (nrow(X_fir_cond_clean) != length(ybar_clean)){
      stop(sprintf("Critical internal error (estimate_hrf_for_condition): Mismatch in length of ybar_clean (%d) and rows of X_fir_cond_clean (%d) for condition: %s",
                 length(ybar_clean), nrow(X_fir_cond_clean), condition_name))
  }
  
  if (sum(abs(X_fir_cond_clean)) < 1e-6 ) { 
      warning(sprintf("No effective stimulus events for condition '%s' fall within non-spike TRs for FIR estimation (in estimate_hrf_for_condition).", condition_name))
      return(current_result)
  }

  min_trs_needed <- max(user_options$hrf_fir_taps * 2, user_options$cv_folds + 1) 
  if (nrow(X_fir_cond_clean) < min_trs_needed) { 
      warning(sprintf("Not enough valid TRs (%d, need ~%d) to estimate FIR for condition: %s (in estimate_hrf_for_condition).", 
                      nrow(X_fir_cond_clean), min_trs_needed, condition_name))
      return(current_result)
  }
  
  cv_results <- cv_fusedlasso(y = ybar_clean, X = X_fir_cond_clean,
                                lambda_grid = user_options$lambda1_grid, 
                                gamma_grid = user_options$lambda2_grid,  
                                k_folds = user_options$cv_folds,
                                block_ids = block_ids_for_cv)
  
  best_lambda1 <- cv_results$best_lambda
  best_lambda2 <- cv_results$best_gamma
  
  message(sprintf("  Condition '%s': Best lambda1 (L1 diff): %.4f, Best lambda2 (L2 coeff): %.4f", 
                  condition_name, best_lambda1, best_lambda2))
  
  X_mat <- X_fir_cond_clean
  p_coeffs_final <- ncol(X_mat)
  D_1d_final <- NULL
  if (p_coeffs_final > 1) {
    D_1d_final <- diff(diag(p_coeffs_final), differences = 1)
  } else if (p_coeffs_final == 1) {
    D_1d_final <- diff(diag(1), differences = 1)
  }
  
  final_fit_path <- tryCatch({
    genlasso::fusedlasso(y = ybar_clean, X = X_mat, D = D_1d_final, gamma = best_lambda2)
  }, error = function(e) {
    warning(sprintf("Final genlasso::fusedlasso fit failed for condition '%s' with gamma=%.4f: %s", 
                    condition_name, best_lambda2, e$message))
    return(NULL)
  })
  
  if (is.null(final_fit_path) || !inherits(final_fit_path, "genlasso")) {
    warning(sprintf("Final genlasso::fusedlasso fit for condition '%s' is NULL or not a genlasso object. Check warnings.", condition_name))
    return(current_result)
  }
  
  fir_coeffs <- tryCatch({
    stats::coef(final_fit_path, lambda = best_lambda1)$beta
  }, error = function(e) {
    warning(sprintf("Failed to extract coefficients for condition '%s' with lambda=%.4f from genlasso fit: %s", 
                    condition_name, best_lambda1, e$message))
    return(NULL)
  })
  
  if (is.null(fir_coeffs)) {
    return(current_result) # Return NULL hrf_estimate as fir_coeffs extraction failed
  }
  
  if (length(fir_coeffs) != user_options$hrf_fir_taps) {
    warning(sprintf("Coefficient length mismatch for condition %s. Expected: %d, Got: %d. Padding/truncating.",
                    condition_name, user_options$hrf_fir_taps, length(fir_coeffs)))
    temp_coeffs <- rep(0, user_options$hrf_fir_taps)
    len_to_copy <- min(length(fir_coeffs), user_options$hrf_fir_taps)
    if (len_to_copy > 0) temp_coeffs[1:len_to_copy] <- fir_coeffs[1:len_to_copy]
    fir_coeffs <- temp_coeffs
  }
  
  current_result$hrf_estimate <- project_cone_heuristic(fir_coeffs)
  current_result$taps <- seq_len(user_options$hrf_fir_taps)
  if (user_options$return_full_model) {
      current_result$glmgen_fit <- final_fit_path
  }
  return(current_result)
} 
