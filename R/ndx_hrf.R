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
#'   - `good_voxel_R2_threshold` (numeric): R-squared threshold for selecting good voxels
#'       (e.g., 0.02-0.06). Use `-Inf` to disable voxel filtering.
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
  validated_spike_TR_mask <- validated_inputs$validated_spike_TR_mask

  # Prepare data for good voxels and clustering
  response_data <- prepare_hrf_response_data(Y_fmri, pass0_residuals, run_idx, 
                                             validated_spike_TR_mask, 
                                             user_options)
  if (is.null(response_data)) {
    warning("HRF estimation cannot proceed: prepare_hrf_response_data returned NULL (e.g. no good voxels or valid TRs).")
    return(NULL) 
  }
  ybar_global_clean <- response_data$ybar_clean # Global ybar for fallback or C=1
  block_ids_for_cv_global <- response_data$block_ids_for_cv
  valid_TRs_for_hrf_estimation_mask <- response_data$valid_TRs_for_hrf_estimation_mask
  good_voxels_indices_in_Y_fmri <- response_data$good_voxels_indices # Assuming prepare_hrf_response_data returns this
  Y_good_voxels_all_trs <- Y_fmri[, good_voxels_indices_in_Y_fmri, drop = FALSE]

  # Perform clustering if method is not "none" and num_clusters > 1
  cluster_info <- NULL
  if (user_options$hrf_cluster_method == "pam" && user_options$num_hrf_clusters > 1 && !is.null(Y_good_voxels_all_trs) && ncol(Y_good_voxels_all_trs) > 0) {
    Y_good_voxels_for_clustering <- Y_good_voxels_all_trs[valid_TRs_for_hrf_estimation_mask, , drop = FALSE]
    if (nrow(Y_good_voxels_for_clustering) > 0) {
      cluster_info <- .perform_hrf_clustering(Y_good_voxels_for_clustering, 
                                            user_options$num_hrf_clusters, 
                                            user_options, verbose = (user_options$verbose_hrf %||% FALSE))
    } else {
        if(user_options$verbose_hrf %||% FALSE) message("No valid TRs in Y_good_voxels_for_clustering, skipping clustering.")
    }
  }
  
  if (is.null(cluster_info) || cluster_info$num_clusters_eff == 1) {
    # Treat as a single cluster (global HRF estimation)
    num_effective_clusters <- 1
    voxel_to_cluster_map <- rep(1L, ncol(Y_good_voxels_all_trs %||% matrix(ncol=0))) # Map good voxels to cluster 1
  } else {
    num_effective_clusters <- cluster_info$num_clusters_eff
    voxel_to_cluster_map <- cluster_info$voxel_cluster_ids # Maps good_voxels_indices to cluster_id
  }

  results_list_collector <- list()
  conditions <- unique(as.character(events$condition))
  
  unique_runs_overall <- unique(run_idx)
  run_lengths_overall <- as.numeric(table(factor(run_idx, levels = unique_runs_overall)))
  sf_overall <- fmrireg::sampling_frame(blocklens = run_lengths_overall, TR = TR)
  # message(sprintf("[ndx_estimate_initial_hrfs] sf_overall$total_samples: %s (class: %s)", 
  #                 as.character(sum(sf_overall$blocklens)), class(sum(sf_overall$blocklens))))
  # if (is.null(sum(sf_overall$blocklens))) {message("sf_overall$total_samples IS NULL - THIS IS THE PROBLEM"); print(sf_overall)}

  for (cond_name_iter in conditions) {
    if (user_options$verbose_hrf %||% FALSE) message(paste("  Processing HRF for condition:", cond_name_iter))
    events_for_this_condition <- events[events$condition == cond_name_iter, , drop = FALSE]

    for (cl_idx in 1:num_effective_clusters) {
      ybar_for_this_hrf <- NULL
      block_ids_for_this_hrf_cv <- NULL

      if (num_effective_clusters == 1) {
        if (user_options$verbose_hrf %||% FALSE) message(sprintf("    Condition %s: Using global ybar for estimation.", cond_name_iter))
        ybar_for_this_hrf <- ybar_global_clean
        block_ids_for_this_hrf_cv <- block_ids_for_cv_global
      } else {
        vox_in_this_cluster_mask <- (voxel_to_cluster_map == cl_idx)
        if (sum(vox_in_this_cluster_mask) == 0) {
          if (user_options$verbose_hrf %||% FALSE) message(sprintf("    Condition %s, Cluster %d: No voxels. Skipping.", cond_name_iter, cl_idx))
          next
        }
        Y_cluster_data <- Y_good_voxels_all_trs[, vox_in_this_cluster_mask, drop = FALSE]
        ybar_cluster_all_trs <- matrixStats::rowMedians(Y_cluster_data, na.rm = TRUE)
        ybar_for_this_hrf <- ybar_cluster_all_trs[valid_TRs_for_hrf_estimation_mask]
        block_ids_for_this_hrf_cv <- run_idx[valid_TRs_for_hrf_estimation_mask]
        if (user_options$verbose_hrf %||% FALSE) message(sprintf("    Condition %s, Cluster %d: Using ybar from %d voxels.", 
                                                               cond_name_iter, cl_idx, sum(vox_in_this_cluster_mask)))
      }

      if (length(ybar_for_this_hrf) == 0) {
          if (user_options$verbose_hrf %||% FALSE) message("      ybar_for_this_hrf is empty. Skipping this HRF estimation.")
          next
      }

      condition_cluster_hrf_result <- estimate_hrf_for_condition(
        condition_name = cond_name_iter,
        events_for_condition = events_for_this_condition,
        ybar_clean = ybar_for_this_hrf,
        block_ids_for_cv = block_ids_for_this_hrf_cv,
        overall_sampling_frame = sf_overall,
        valid_TRs_mask = valid_TRs_for_hrf_estimation_mask, 
        TR = TR,
        user_options = user_options # Pass full user_options (contains opts_hrf implicitly)
      )
      
      # Store result with cluster_id if multiple clusters, otherwise cluster_id = 1
      # The estimate_hrf_for_condition already returns a list with $condition, $hrf_estimate etc.
      # We need to add cluster_id to it before collecting.
      if (!is.null(condition_cluster_hrf_result)) {
          condition_cluster_hrf_result$cluster_id <- cl_idx
          results_list_collector[[paste(cond_name_iter, cl_idx, sep="_")]] <- condition_cluster_hrf_result
      }
    }
  }
  
  if (length(results_list_collector) == 0) {
    warning("HRF estimation: No results collected for any condition/cluster.")
    return(NULL)
  }
  
  valid_results <- Filter(Negate(is.null), results_list_collector)
  if(length(valid_results) == 0) {
      warning("No valid HRF results were generated for any condition/cluster.")
      return(NULL)
  }

  # Modify output tibble structure for cluster information
  conditions_out   <- sapply(valid_results, function(x) x$condition)
  cluster_ids_out  <- sapply(valid_results, function(x) x$cluster_id)
  hrf_estimates_out <- lapply(valid_results, function(x) x$hrf_estimate)
  taps_out         <- lapply(valid_results, function(x) x$taps)
  
  output_tbl <- tibble::tibble(
      condition    = conditions_out,
      cluster_id   = cluster_ids_out, 
      hrf_estimate = hrf_estimates_out,
      taps         = taps_out
  )
  
  if (user_options$return_full_model) {
      glmgen_fits_out <- lapply(valid_results, function(x) x$glmgen_fit %||% list(NULL)) # Ensure list if fit is NULL
      if (nrow(output_tbl) > 0) {
        output_tbl$glmgen_fit <- glmgen_fits_out
      } 
  }
  
  if (nrow(output_tbl) == 0 && length(conditions) > 0) {
      warning("HRF estimation resulted in an empty table despite conditions being present. Check warnings.")
      return(NULL)
  }

  # Add attribute for how many clusters were effectively used/estimated
  attr(output_tbl, "num_effective_clusters") <- num_effective_clusters
  if (!is.null(cluster_info) && !is.null(cluster_info$medoid_indices)) {
      # Map medoid_indices (which are relative to good_voxels_indices) back to original Y_fmri voxel indices
      attr(output_tbl, "medoid_voxel_original_indices") <- good_voxels_indices_in_Y_fmri[cluster_info$medoid_indices]
  }
  if (!is.null(voxel_to_cluster_map)) {
      # This map is for good_voxels_indices_in_Y_fmri -> cluster_id
      # To make it map original voxel_id -> cluster_id:
      full_voxel_to_cluster_map <- rep(NA_integer_, ncol(Y_fmri))
      if(length(good_voxels_indices_in_Y_fmri) == length(voxel_to_cluster_map)) {
          full_voxel_to_cluster_map[good_voxels_indices_in_Y_fmri] <- voxel_to_cluster_map
      }
      attr(output_tbl, "voxel_to_cluster_assignment_full") <- full_voxel_to_cluster_map
  }

  return(output_tbl)
}

#' Generate FIR Design Matrix for a Specific Condition
#' @keywords internal
get_fir_design_matrix_for_condition <- function(condition_name, events_df,
                                                sampling_frame, fir_taps, TR) {
  
  total_timepoints <- sum(sampling_frame$blocklens)
  
  # --- BEGIN DEBUG MESSAGES ---
  message(sprintf("[get_fir_design_matrix_for_condition] Cond: %s", condition_name))
  message(sprintf("  Input fir_taps: %s (class: %s), Input TR: %s (class: %s)", 
                  as.character(fir_taps), class(fir_taps), as.character(TR), class(TR)))
  message(sprintf("  sampling_frame$total_samples (total_timepoints): %s (class: %s)", 
                  as.character(total_timepoints), class(total_timepoints)))
  if (!is.numeric(total_timepoints) || length(total_timepoints) != 1 || !is.finite(total_timepoints) || total_timepoints < 0) {
      message("  WARNING: total_timepoints is not a single positive finite number!")
  }
  if (!is.numeric(fir_taps) || length(fir_taps) != 1 || !is.finite(fir_taps) || fir_taps < 0) {
      message("  WARNING: fir_taps is not a single positive finite number!")
  }
  # --- END DEBUG MESSAGES ---

  if (nrow(events_df) == 0 || fir_taps == 0) {
    # If no events or no FIR taps requested, return a zero matrix of correct dimensions
    # This is important for consistency if other conditions *do* produce regressors.
    return(matrix(0, nrow = total_timepoints, ncol = fir_taps))
  }
  
  # Manual FIR basis matrix construction (similar to .ndx_generate_task_regressors)
  X_fir_basis_cond <- matrix(0, nrow = total_timepoints, ncol = fir_taps)
  fir_colnames <- paste0(make.names(condition_name), "_FIRbasis", 0:(fir_taps-1)) # Unique basis names
  colnames(X_fir_basis_cond) <- fir_colnames

  current_run_offset <- 0 # To handle timepoints across concatenated runs
  for (run_idx_val in 1:length(sampling_frame$blocklens)) {
    run_length_tps <- sampling_frame$blocklens[run_idx_val]
    run_start_tp_global <- current_run_offset + 1
    run_end_tp_global <- current_run_offset + run_length_tps
    
    run_events <- events_df[events_df$blockids == run_idx_val, , drop = FALSE]

    if (nrow(run_events) > 0) {
      for (ev_idx in 1:nrow(run_events)) {
        event_onset_sec <- run_events$onsets[ev_idx]
        event_onset_tr_global_0idx <- floor(event_onset_sec / TR) + (run_start_tp_global - 1)

        for (tap_idx in 0:(fir_taps - 1)) { # 0-indexed tap
          target_timepoint_0idx <- event_onset_tr_global_0idx + tap_idx
          target_timepoint_1idx <- target_timepoint_0idx + 1

          if (target_timepoint_1idx >= run_start_tp_global && 
              target_timepoint_1idx <= run_end_tp_global && 
              target_timepoint_1idx <= total_timepoints) {
            X_fir_basis_cond[target_timepoint_1idx, tap_idx + 1] <- 1
          }
        }
      }
    }
    current_run_offset <- current_run_offset + run_length_tps
  }
  
  # The original function had checks for ncol(X_fir) != fir_taps
  # With manual construction, this should always match if fir_taps > 0.
  # If fir_taps was 0, we returned a 0-col matrix earlier.
  if (ncol(X_fir_basis_cond) != fir_taps) {
      # This case should ideally not be reached if fir_taps > 0 due to pre-allocation
      warning(sprintf("Manual FIR design for %s has %d columns, expected %d. This is unexpected.", 
                      condition_name, ncol(X_fir_basis_cond), fir_taps))
      # Fallback to a zero matrix of correct dimensions if something went very wrong
      return(matrix(0, nrow = total_timepoints, ncol = fir_taps, dimnames=list(NULL, fir_colnames))) 
  }
  
  return(X_fir_basis_cond)
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
    return_full_model = FALSE,
    hrf_cluster_method = "none", 
    num_hrf_clusters = 1L,
    hrf_min_events_for_fir = 6L, 
    hrf_low_event_threshold = 12L,
    hrf_target_event_count_for_lambda_scaling = 20L,
    hrf_use_canonical_fallback_for_ultra_sparse = FALSE
  )
  user_options <- utils::modifyList(default_opts, user_options)

  # Check for presence of all expected options now includes new ones
  expected_opts_names <- names(default_opts)
  if(!all(expected_opts_names %in% names(user_options))) {
    stop(paste("Missing one or more user_options for HRF estimation (after merging defaults):", 
               paste(expected_opts_names[!expected_opts_names %in% names(user_options)], collapse=", ")))
  }

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
  
  # Validate specific user_options types and values (add new ones)
  if(!is.numeric(user_options$hrf_fir_taps) || user_options$hrf_fir_taps <=0 || round(user_options$hrf_fir_taps) != user_options$hrf_fir_taps) stop("hrf_fir_taps must be a positive integer.")
  user_options$hrf_fir_taps <- as.integer(user_options$hrf_fir_taps)

  if(!is.numeric(user_options$hrf_fir_span_seconds) || user_options$hrf_fir_span_seconds <=0) stop("hrf_fir_span_seconds must be a positive number.")
  
  thr <- user_options$good_voxel_R2_threshold
  if(!is.numeric(thr) || length(thr) != 1 || is.na(thr) ||
     (is.finite(thr) && (thr < -1 || thr > 1)) || # Allow -Inf but not other negative numbers for R2 unless it means something
     (is.infinite(thr) && thr > 0)) {
    stop("good_voxel_R2_threshold must be -Inf or a finite number between 0 and 1 (or slightly <0 if that has meaning).")
  }
  
  if(!is.numeric(user_options$cv_folds) || user_options$cv_folds <=0 || round(user_options$cv_folds) != user_options$cv_folds) stop("cv_folds must be a positive integer.")
  user_options$cv_folds <- as.integer(user_options$cv_folds)

  if(!is.numeric(user_options$lambda1_grid) || length(user_options$lambda1_grid)==0) stop("lambda1_grid must be a non-empty numeric vector.")
  if(!is.numeric(user_options$lambda2_grid) || length(user_options$lambda2_grid)==0) stop("lambda2_grid must be a non-empty numeric vector.")
  
  if(!is.numeric(user_options$hrf_min_good_voxels) || user_options$hrf_min_good_voxels <0 || round(user_options$hrf_min_good_voxels) != user_options$hrf_min_good_voxels) stop("hrf_min_good_voxels must be a non-negative integer.")
  user_options$hrf_min_good_voxels <- as.integer(user_options$hrf_min_good_voxels)
  
  if(!is.logical(user_options$return_full_model) || length(user_options$return_full_model)!=1) stop("return_full_model must be a single logical value.")

  if(!is.character(user_options$hrf_cluster_method) || length(user_options$hrf_cluster_method)!=1 || !user_options$hrf_cluster_method %in% c("none", "pam")) stop("hrf_cluster_method must be 'none' or 'pam'.")
  
  if(!is.numeric(user_options$num_hrf_clusters) || user_options$num_hrf_clusters < 1 || round(user_options$num_hrf_clusters) != user_options$num_hrf_clusters) stop("num_hrf_clusters must be a positive integer >= 1.")
  user_options$num_hrf_clusters <- as.integer(user_options$num_hrf_clusters)
  if (user_options$hrf_cluster_method == "none" && user_options$num_hrf_clusters != 1) {
      warning("hrf_cluster_method is 'none' but num_hrf_clusters is not 1. Setting num_hrf_clusters to 1.")
      user_options$num_hrf_clusters <- 1L
  }

  if(!is.numeric(user_options$hrf_min_events_for_fir) || user_options$hrf_min_events_for_fir < 0 || round(user_options$hrf_min_events_for_fir) != user_options$hrf_min_events_for_fir) stop("hrf_min_events_for_fir must be a non-negative integer.")
  user_options$hrf_min_events_for_fir <- as.integer(user_options$hrf_min_events_for_fir)

  if(!is.numeric(user_options$hrf_low_event_threshold) || user_options$hrf_low_event_threshold < 0 || round(user_options$hrf_low_event_threshold) != user_options$hrf_low_event_threshold) stop("hrf_low_event_threshold must be a non-negative integer.")
  user_options$hrf_low_event_threshold <- as.integer(user_options$hrf_low_event_threshold)
  if(user_options$hrf_low_event_threshold < user_options$hrf_min_events_for_fir) {
      warning("hrf_low_event_threshold should not be less than hrf_min_events_for_fir. Adjusting hrf_low_event_threshold.")
      user_options$hrf_low_event_threshold <- user_options$hrf_min_events_for_fir
  }

  if(!is.numeric(user_options$hrf_target_event_count_for_lambda_scaling) || user_options$hrf_target_event_count_for_lambda_scaling <=0 || round(user_options$hrf_target_event_count_for_lambda_scaling) != user_options$hrf_target_event_count_for_lambda_scaling) stop("hrf_target_event_count_for_lambda_scaling must be a positive integer.")
  user_options$hrf_target_event_count_for_lambda_scaling <- as.integer(user_options$hrf_target_event_count_for_lambda_scaling)
  
  if(!is.logical(user_options$hrf_use_canonical_fallback_for_ultra_sparse) || length(user_options$hrf_use_canonical_fallback_for_ultra_sparse)!=1) stop("hrf_use_canonical_fallback_for_ultra_sparse must be a single logical value.")

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
  
  ybar_all_trs <- matrixStats::rowMedians(Y_good_voxels, na.rm = TRUE) 
  ybar_clean <- ybar_all_trs[valid_TRs_for_hrf_estimation_mask]
  block_ids_for_cv <- run_idx[valid_TRs_for_hrf_estimation_mask]
  
  return(list(
    ybar_clean = ybar_clean,
    block_ids_for_cv = block_ids_for_cv,
    valid_TRs_for_hrf_estimation_mask = valid_TRs_for_hrf_estimation_mask,
    n_good_voxels = n_good_voxels,
    good_voxels_indices = good_voxels_idx
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
  
  opts_hrf <- user_options # If user_options is already opts_hrf, then this is fine.
                          # If user_options is the global one, then opts_hrf <- user_options$opts_hrf
                          # Assuming it's already the HRF-specific sub-list from the caller.

  current_result <- list(condition = condition_name, hrf_estimate = NULL, taps = NULL, glmgen_fit = NULL)
  # if (opts_hrf$return_full_model) current_result$glmgen_fit <- NULL # Already init with NULL

  if (nrow(events_for_condition) == 0) {
      # This warning might be redundant if ndx_estimate_initial_hrfs already filters conditions with no events.
      # However, events_for_condition here is specific to current cluster's ybar context for valid TRs.
      warning(paste("No events found for condition ", condition_name, " to estimate HRF.")) 
      return(current_result)
  }
  
  # Calculate n_events_q for the current context (valid TRs for ybar_clean)
  # This needs the full events table and the valid_TRs_mask to count events whose time window falls into valid TRs.
  # For simplicity now, we use nrow(events_for_condition) as a proxy, assuming all these events contribute.
  # A more precise n_events_q would consider event durations and overlap with valid_TRs_mask.
  n_events_q <- nrow(events_for_condition) 
  lambda_scale_factor <- 1.0
  effective_hrf_fir_taps <- opts_hrf$hrf_fir_taps

  if (n_events_q < (opts_hrf$hrf_min_events_for_fir %||% 6L)) {
    if (opts_hrf$hrf_use_canonical_fallback_for_ultra_sparse %||% FALSE) {
      warning(sprintf("Condition '%s' has only %d events. Using canonical HRF (not yet implemented, returning zeros).", condition_name, n_events_q))
      # Placeholder: return a zero HRF or a pre-defined canonical scaled by some basic fit
      current_result$hrf_estimate <- rep(0, opts_hrf$hrf_fir_taps)
      current_result$taps <- seq_len(opts_hrf$hrf_fir_taps)
      return(current_result)
    } else {
      warning(sprintf("Condition '%s' has only %d events (min_events_for_fir=%d). Returning zero/damped FIR.", 
                      condition_name, n_events_q, (opts_hrf$hrf_min_events_for_fir %||% 6L)))
      current_result$hrf_estimate <- rep(0, opts_hrf$hrf_fir_taps) # Zero HRF
      current_result$taps <- seq_len(opts_hrf$hrf_fir_taps)
      return(current_result)
    }
  } else if (n_events_q < (opts_hrf$hrf_low_event_threshold %||% 12L)) {
    target_events <- (opts_hrf$hrf_target_event_count_for_lambda_scaling %||% 20L)
    if (n_events_q > 0) { # Avoid division by zero if n_events_q is somehow zero here despite previous check
        lambda_scale_factor <- sqrt(target_events / n_events_q)
    }
    if (opts_hrf$verbose_hrf %||% FALSE) message(sprintf("  Condition '%s' has %d events. Applying lambda scaling factor: %.2f", 
                                 condition_name, n_events_q, lambda_scale_factor))
  }

  X_fir_cond_all_trs <- tryCatch({
       get_fir_design_matrix_for_condition(
          condition_name = condition_name,
          events_df = events_for_condition, 
          sampling_frame = overall_sampling_frame, 
          fir_taps = effective_hrf_fir_taps, # Use effective_hrf_fir_taps
          TR = TR 
      )
  }, error = function(e) {
      warning(paste("Error generating FIR design for condition", condition_name, "in estimate_hrf_for_condition:", e$message))
      return(NULL)
  })

  if (is.null(X_fir_cond_all_trs) || ncol(X_fir_cond_all_trs) == 0) {
    warning(paste("FIR design matrix for condition", condition_name, "is empty or NULL."))
    return(current_result)
  }
  
  X_fir_cond_clean <- X_fir_cond_all_trs[valid_TRs_mask, , drop = FALSE]
  
  if (nrow(X_fir_cond_clean) != length(ybar_clean)){
      stop(sprintf("Critical internal error: Mismatch in length of ybar_clean (%d) and rows of X_fir_cond_clean (%d) for condition: %s",
                 length(ybar_clean), nrow(X_fir_cond_clean), condition_name))
  }
  
  if (sum(abs(X_fir_cond_clean)) < 1e-6 ) { 
      warning(sprintf("No effective stimulus events for condition '%s' fall within non-spike TRs for FIR estimation.", condition_name))
      return(current_result)
  }

  min_trs_needed <- max(effective_hrf_fir_taps * 2, (opts_hrf$cv_folds %||% 2L) + 1) 
  if (nrow(X_fir_cond_clean) < min_trs_needed) { 
      warning(sprintf("Not enough valid TRs (%d, need ~%d) for FIR for condition: %s.", 
                      nrow(X_fir_cond_clean), min_trs_needed, condition_name))
      return(current_result)
  }
  
  # Apply scaling to lambda grids for CV if lambda_scale_factor > 1
  # Note: genlasso takes lambda (for L1 on diffs) and gamma (for L2 on coeffs)
  # The feedback suggested scaling both lambda1 and lambda2 (gamma).
  cv_lambda1_grid <- (opts_hrf$lambda1_grid %||% 1) * lambda_scale_factor
  cv_lambda2_grid <- (opts_hrf$lambda2_grid %||% 0.1) * lambda_scale_factor # Scale gamma as well

  cv_results <- cv_fusedlasso(y = ybar_clean, X = X_fir_cond_clean,
                                lambda_grid = cv_lambda1_grid, 
                                gamma_grid = cv_lambda2_grid,  
                                k_folds = (opts_hrf$cv_folds %||% 2L),
                                block_ids = block_ids_for_cv)
  
  # Use the best lambdas directly as returned by cv_fusedlasso (which were sought on potentially scaled grids)
  best_lambda1 <- cv_results$best_lambda 
  best_lambda2 <- cv_results$best_gamma  
  
  if (opts_hrf$verbose_hrf %||% FALSE) {
    message(sprintf("  Condition '%s': Best lambda1 (L1 diff): %.4f, Best lambda2 (L2 coeff): %.4f (from scaled grids if applicable, scale factor: %.2f)", 
                  condition_name, best_lambda1, best_lambda2, lambda_scale_factor))
  }
  
  X_mat <- X_fir_cond_clean
  p_coeffs_final <- ncol(X_mat)
  D_1d_final <- NULL
  if (p_coeffs_final > 1) {
    D_1d_final <- diff(diag(p_coeffs_final), differences = 1)
  } else if (p_coeffs_final == 1) {
    # For genlasso, if D is NULL or not matrix, it implies no fusion. For single coeff, this is fine.
    D_1d_final <- matrix(0, nrow=0, ncol=1) # genlasso needs D to be matrix, even if 0-row for single coeff
  }
  
  final_fit_path <- tryCatch({
    genlasso::fusedlasso(y = ybar_clean, X = X_mat, D = D_1d_final, gamma = best_lambda2)
  }, error = function(e) {
    warning(sprintf("Final genlasso::fusedlasso fit failed for condition '%s' with gamma=%.4f: %s", 
                    condition_name, best_lambda2, e$message))
    return(NULL)
  })
  
  if (is.null(final_fit_path) || !inherits(final_fit_path, "genlasso")) {
    warning(sprintf("Final genlasso::fusedlasso fit for condition '%s' is NULL or not a genlasso object.", condition_name))
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
    return(current_result) 
  }
  
  # Ensure fir_coeffs length matches effective_hrf_fir_taps
  if (length(fir_coeffs) != effective_hrf_fir_taps) {
    if (opts_hrf$verbose_hrf %||% FALSE) {
        message(sprintf("    Condition %s: Coeff length %d from fusedlasso, expected %d. Padding/truncating.",
                        condition_name, length(fir_coeffs), effective_hrf_fir_taps)) 
    }
    temp_coeffs <- rep(0, effective_hrf_fir_taps)
    len_to_copy <- min(length(fir_coeffs), effective_hrf_fir_taps)
    if (len_to_copy > 0) temp_coeffs[1:len_to_copy] <- fir_coeffs[1:len_to_copy]
    fir_coeffs <- temp_coeffs
  }
  
  current_result$hrf_estimate <- project_cone_heuristic(fir_coeffs)
  current_result$taps <- seq_len(effective_hrf_fir_taps)
  if (opts_hrf$return_full_model) {
      current_result$glmgen_fit <- final_fit_path
  }
  return(current_result)
}

#' Perform K-Medoids Clustering on Voxel Time Courses
#' @keywords internal
.perform_hrf_clustering <- function(Y_for_clustering, num_clusters, user_options, verbose) {
  # Y_for_clustering is timepoints x good_voxels_for_clustering
  # We cluster voxels, so pam needs voxels x timepoints, or dist on voxels
  if (verbose) message(sprintf("  Performing K-Medoids clustering for %d clusters on %d voxels.", num_clusters, ncol(Y_for_clustering)))
  
  if (num_clusters == 1 || ncol(Y_for_clustering) <= num_clusters || ncol(Y_for_clustering) < (user_options$hrf_min_good_voxels %||% 10) ) {
    if (verbose && num_clusters > 1) message("    Number of good voxels too small for requested clusters, or num_clusters=1. Assigning all to cluster 1.")
    return(list(voxel_cluster_ids = rep(1L, ncol(Y_for_clustering)), 
                num_clusters_eff = 1L, 
                medoid_indices = if(ncol(Y_for_clustering)>0) 1L else integer(0) )) # Default to first voxel as medoid if any
  }
  
  # cluster::pam expects samples in rows, features in columns. 
  # To cluster voxels (features here), we can transpose Y_for_clustering.
  # Or, supply a distance matrix. Transposing is simpler if data fits memory.
  # Consider if scaling/normalizing voxels before clustering is needed (e.g. z-score each voxel's timecourse)
  pam_fit <- NULL
  tryCatch({
    # Note: pam can be slow for large N (voxels) and large P (timepoints per voxel)
    # For very large number of voxels, might need to subsample or use clara.
    # Data for pam should be n_observations (voxels) x n_variables (timepoints)
    pam_fit <- cluster::pam(t(Y_for_clustering), k = num_clusters, diss = FALSE, metric="euclidean", stand = FALSE)
  }, error = function(e) {
    warning(sprintf("K-Medoids clustering (pam) failed: %s. Assigning all to cluster 1.", e$message))
    pam_fit <<- NULL # Ensure it's NULL in outer scope on error
  })
  
  if (is.null(pam_fit)) {
    return(list(voxel_cluster_ids = rep(1L, ncol(Y_for_clustering)), 
                num_clusters_eff = 1L, 
                medoid_indices = if(ncol(Y_for_clustering)>0) 1L else integer(0) ))
  }
  
  return(list(voxel_cluster_ids = pam_fit$clustering, 
              num_clusters_eff = length(unique(pam_fit$clustering)), 
              medoid_indices = pam_fit$id.med # Indices of medoids *within the input to pam* (i.e., among good_voxels_idx)
              ))
} 
