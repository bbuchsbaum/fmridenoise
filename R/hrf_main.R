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
#'                                  motion_params_example, run_idx_example, TR,
#'                                  poly_degree = 1L)
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
      raw_cluster_info <- .perform_hrf_clustering(Y_good_voxels_for_clustering,
                                            user_options$num_hrf_clusters,
                                            user_options, verbose = (user_options$verbose_hrf %||% FALSE))
      cluster_info <- .auto_adapt_hrf_clusters(
                         Y_for_clustering = Y_good_voxels_for_clustering,
                         cluster_ids = raw_cluster_info$voxel_cluster_ids,
                         merge_corr_thresh = user_options$hrf_cluster_merge_corr_thresh,
                         min_cluster_size = user_options$hrf_cluster_min_size,
                         verbose = (user_options$verbose_hrf %||% FALSE))
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
