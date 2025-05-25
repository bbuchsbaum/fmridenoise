#' Assemble Default ND-X User Options
#'
#' Returns a list containing the default option structure used by
#' \code{NDX_Process_Subject}. The values mirror the defaults documented
#' in \code{\link{ndx_user_options}}.
#'
#' @return Named list suitable for the \code{user_options} argument of
#'   \code{NDX_Process_Subject}.
#' @export
ndx_default_user_options <- function() {
  list(
    max_passes = 3L,
    min_des_gain_convergence = 0.005,
    min_rho_noise_projection_convergence = 0.01,
    task_regressor_names_for_extraction = character(0),
    verbose = TRUE,
    opts_pass0 = list(poly_degree = 1),
    opts_hrf = list(
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
    ),
    opts_rpca = list(),
    opts_spectral = list(),
    opts_whitening = list(),
    opts_ridge = list(
      lambda_ridge = 1.0,
      anisotropic_ridge_enable = TRUE,
      lambda_signal = 0.1,
      lambda_noise_rpca = 1.0,
      lambda_noise_spectral = 1.0,
      lambda_noise_other = 1.0,
      ridge_gcv_folds = 5L,
      ridge_update_lambda_aggressiveness = FALSE,
      lambda_noise_gdlite = 1.0,
      lambda_noise_ndx_unique = 1.0
    ),
    opts_annihilation = list(
      annihilation_enable_mode = FALSE,
      annihilation_gdlite_poly_degree = 1,
      annihilation_gdlite_k_max = 30L,
      annihilation_gdlite_r2_thresh_noise_pool = 0.05,
      annihilation_gdlite_tsnr_thresh_noise_pool = 30,
      annihilation_gdlite_r2_thresh_good_voxels = 0.05
    )
  )
}
