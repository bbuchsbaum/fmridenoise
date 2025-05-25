#' Helper Utilities for ND-X Workflow
#'
#' These functions assist `NDX_Process_Subject` by providing
#' input validation and option handling utilities.
#'
#' @name ndx_workflow_utils
NULL

#' Validate inputs for `NDX_Process_Subject`
#'
#' @param Y_fmri Numeric matrix of fMRI data.
#' @param events Data frame of events.
#' @param motion_params Numeric matrix of motion parameters.
#' @param run_idx Numeric vector indicating run membership.
#' @param TR Numeric repetition time.
#' @return Invisible TRUE if validation passes, otherwise an error is thrown.
#' @export
ndx_validate_process_subject_inputs <- function(Y_fmri, events, motion_params,
                                                run_idx, TR) {
  if (!is.matrix(Y_fmri) || !is.numeric(Y_fmri)) {
    stop("Y_fmri must be a numeric matrix (timepoints x voxels).")
  }
  n_timepoints <- nrow(Y_fmri)
  if (n_timepoints == 0) stop("Y_fmri has zero timepoints.")
  if (ncol(Y_fmri) == 0) stop("Y_fmri has zero voxels.")

  if (!is.matrix(motion_params) || !is.numeric(motion_params)) {
    stop("motion_params must be a numeric matrix.")
  }
  if (nrow(motion_params) != n_timepoints) {
    stop("Number of rows in motion_params must match Y_fmri.")
  }

  if (!is.data.frame(events)) stop("events must be a data frame.")
  required_event_cols <- c("onsets", "durations", "condition", "blockids")
  if (!all(required_event_cols %in% names(events))) {
    stop(paste("events data frame must contain columns:",
               paste(required_event_cols, collapse = ", ")))
  }

  if (!is.numeric(run_idx) || length(run_idx) != n_timepoints) {
    stop("run_idx must be a numeric vector with length matching nrow(Y_fmri).")
  }
  if (!is.numeric(TR) || length(TR) != 1 || TR <= 0) {
    stop("TR must be a single positive number.")
  }
  invisible(TRUE)
}

#' Merge user options with defaults
#'
#' @param user_options List of user supplied options.
#' @return Named list of options where unspecified values are filled with defaults.
#' @export
ndx_prepare_workflow_options <- function(user_options = list()) {
  defaults <- ndx_default_user_options()
  utils::modifyList(defaults, user_options)
}

#' Run the optional GLMdenoise-Lite setup for Annihilation mode
#'
#' @param Y_fmri Numeric matrix of fMRI data.
#' @param events Data frame of events.
#' @param motion_params Motion parameters matrix.
#' @param run_idx Run index vector.
#' @param TR Repetition time in seconds.
#' @param opts_annihilation List of annihilation options.
#' @param verbose Logical verbosity flag.
#' @return List with `selected_pcs` and full `gdlite_results` if run.
#' @export
ndx_run_annihilation_setup <- function(Y_fmri, events, motion_params, run_idx, TR,
                                       opts_annihilation, verbose = TRUE) {
  result <- list(selected_pcs = NULL, gdlite_results = NULL)
  if (!isTRUE(opts_annihilation$annihilation_enable_mode)) {
    return(result)
  }

  if (verbose) message("Annihilation Mode enabled: Running GLMdenoise-Lite...")
  original <- getOption("fmridenoise.verbose_gdlite")
  options(fmridenoise.verbose_gdlite = verbose)
  gdlite_initial_results <- ndx_run_gdlite(
    Y_fmri = Y_fmri,
    events = events,
    run_idx = run_idx,
    TR = TR,
    motion_params = motion_params,
    poly_degree = opts_annihilation$annihilation_gdlite_poly_degree,
    k_max = opts_annihilation$annihilation_gdlite_k_max,
    r2_thresh_noise_pool = opts_annihilation$annihilation_gdlite_r2_thresh_noise_pool,
    tsnr_thresh_noise_pool = opts_annihilation$annihilation_gdlite_tsnr_thresh_noise_pool,
    r2_thresh_good_voxels = opts_annihilation$annihilation_gdlite_r2_thresh_good_voxels,
    perform_final_glm = FALSE
  )
  options(fmridenoise.verbose_gdlite = original)

  if (!is.null(gdlite_initial_results) && !is.null(gdlite_initial_results$selected_pcs) &&
      ncol(gdlite_initial_results$selected_pcs) > 0) {
    result$selected_pcs <- gdlite_initial_results$selected_pcs
    if (verbose) message(sprintf("  GLMdenoise-Lite selected %d PCs.",
                                 ncol(result$selected_pcs)))
  } else if (verbose) {
    message("  GLMdenoise-Lite did not select any PCs or failed.")
  }
  result$gdlite_results <- gdlite_initial_results
  result
}
