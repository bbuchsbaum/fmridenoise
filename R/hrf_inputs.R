#' @keywords internal
# NO @export for validate_hrf_inputs
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
    hrf_cluster_merge_corr_thresh = 0.95,
    hrf_cluster_min_size = 5L,
    hrf_min_events_for_fir = 6L,
    hrf_low_event_threshold = 12L,
    hrf_target_event_count_for_lambda_scaling = 20L,
    hrf_use_canonical_fallback_for_ultra_sparse = FALSE,
    hrf_cone_nonneg = TRUE,
    hrf_cone_unimodal = TRUE,
    hrf_cone_normalize_area = TRUE,
    verbose_hrf = FALSE 
  )
  user_options <- utils::modifyList(default_opts, user_options)

  # Check for presence of all expected options (now includes new ones)
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

  blocks_events <- sort(unique(events$blockids))
  blocks_runs <- sort(unique(run_idx))
  if (!identical(blocks_events, blocks_runs)) {
    stop(sprintf(
      "events$blockids %s do not match run_idx %s",
      paste(blocks_events, collapse = ","),
      paste(blocks_runs, collapse = ",")
    ))
  }
  events <- events[order(factor(events$blockids, levels = unique(run_idx))), , drop = FALSE]

  validated_spike_TR_mask <- spike_TR_mask
  if (is.null(validated_spike_TR_mask)) {
    validated_spike_TR_mask <- rep(FALSE, n_timepoints)
  }
  if (length(validated_spike_TR_mask) != n_timepoints || !is.logical(validated_spike_TR_mask)) {
      stop("spike_TR_mask must be a logical vector with length matching n_timepoints in Y_fmri.")
  }
  
  # Validate specific user_options types and values
  if(!is.numeric(user_options$hrf_fir_taps) || user_options$hrf_fir_taps <=0 || round(user_options$hrf_fir_taps) != user_options$hrf_fir_taps) stop("hrf_fir_taps must be a positive integer.")
  user_options$hrf_fir_taps <- as.integer(user_options$hrf_fir_taps)

  if(!is.numeric(user_options$hrf_fir_span_seconds) || user_options$hrf_fir_span_seconds <=0) stop("hrf_fir_span_seconds must be a positive number.")
  
  thr <- user_options$good_voxel_R2_threshold
  if(!is.numeric(thr) || length(thr) != 1 || is.na(thr) ||
     (is.finite(thr) && (thr < -1 || thr > 1)) || 
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

  if(!is.numeric(user_options$hrf_cluster_merge_corr_thresh) || length(user_options$hrf_cluster_merge_corr_thresh) != 1 ||
     user_options$hrf_cluster_merge_corr_thresh <= -1 || user_options$hrf_cluster_merge_corr_thresh > 1) {
      stop("hrf_cluster_merge_corr_thresh must be a single numeric value in (-1,1].")
  }
  if(!is.numeric(user_options$hrf_cluster_min_size) || user_options$hrf_cluster_min_size < 1 || round(user_options$hrf_cluster_min_size) != user_options$hrf_cluster_min_size) {
      stop("hrf_cluster_min_size must be a positive integer.")
  }
  user_options$hrf_cluster_min_size <- as.integer(user_options$hrf_cluster_min_size)

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
  # Validations for new cone projection options
  if(!is.logical(user_options$hrf_cone_nonneg) || length(user_options$hrf_cone_nonneg)!=1) stop("hrf_cone_nonneg must be a single logical value.")
  if(!is.logical(user_options$hrf_cone_unimodal) || length(user_options$hrf_cone_unimodal)!=1) stop("hrf_cone_unimodal must be a single logical value.")
  if(!is.logical(user_options$hrf_cone_normalize_area) || length(user_options$hrf_cone_normalize_area)!=1) stop("hrf_cone_normalize_area must be a single logical value.")
  if(!is.logical(user_options$verbose_hrf) || length(user_options$verbose_hrf)!=1) stop("verbose_hrf must be a single logical value.")

  return(list(user_options = user_options, 
              n_timepoints = n_timepoints, 
              validated_spike_TR_mask = validated_spike_TR_mask))
}
