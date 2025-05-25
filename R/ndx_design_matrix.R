#' Build a Full Design Matrix from Modular Components
#'
#' This is a refactored, cleaner implementation of `ndx_build_design_matrix` that
#' delegates the heavy lifting to a set of small internal helper functions.  The
#' public interface (arguments and return value) is identical to the original
#' monolithic version, so existing code and tests do not need to change.
#'
#' @param estimated_hrfs Tibble of HRF estimates; see old docs.
#' @param events Data frame of experimental events.
#' @param motion_params Optional matrix of motion parameters.
#' @param rpca_components Optional matrix of RPCA regressors.
#' @param spectral_sines Optional matrix of spectral sine/cosine regressors.
#' @param run_idx Integer vector indicating run membership for each time-point.
#' @param TR Positive repetition time, in seconds.
#' @param poly_degree Degree of Legendre polynomial baseline (per run).
#' @param verbose Logical flag for verbose console output.
#' @param drop_zero_variance Logical flag to drop near-zero variance regressors.
#' @param zero_var_epsilon Numeric tolerance for detecting near-zero variance
#'   regressors. Columns with variance below this value are considered
#'   near-constant. Default is 1e-8.
#'
#' @return Numeric matrix with one column per regressor (or `NULL` on failure).
#' @importFrom fmrireg sampling_frame event_model design_matrix baseline_model
#' @importFrom tibble is_tibble
#' @export ndx_build_design_matrix
ndx_build_design_matrix <- function(estimated_hrfs,
                                    events,
                                    motion_params,
                                    rpca_components,
                                    spectral_sines,
                                    run_idx,
                                    TR,
                                    poly_degree = NULL,
                                    verbose = TRUE,
                                    drop_zero_variance = FALSE,
                                    zero_var_epsilon = 1e-8) {

  # 1. Basic validation ----------------------------------------------------
  if (!(is.numeric(TR) && length(TR) == 1 && TR > 0)) {
    stop("TR must be a single positive number.")
  }

  info <- .ndx_validate_design_inputs(run_idx, motion_params, rpca_components,
                                      spectral_sines, events)
  sf   <- fmrireg::sampling_frame(blocklens = info$run_lengths, TR = TR)

  events_mapped <- events
  if (!is.null(events_mapped)) {
    events_mapped$blockids <- info$run_mapping[as.character(events_mapped$blockids)]
  }

  mapped_run_idx <- info$run_idx_mapped


  # 2. Task regressors -----------------------------------------------------
  task_mat <- .ndx_generate_task_regressors(estimated_hrfs, events_mapped, sf, TR, verbose)

  # 3. Nuisance regressors -------------------------------------------------
  nuisance_list_components <- .ndx_generate_nuisance_regressors(motion_params, rpca_components, spectral_sines, verbose)
  
  # 4. Baseline regressors -------------------------------------------------
  baseline_mat <- .ndx_generate_baseline_regressors(mapped_run_idx, sf, poly_degree, verbose)

  # 5. Combine -------------------------------------------------------------
  regressor_list <- c(list(task = task_mat), nuisance_list_components, list(baseline = baseline_mat))
  X_full <- .ndx_combine_regressors(regressor_list, info$total_tp,
                                    verbose, drop_zero_variance,
                                    zero_var_epsilon)

  if (verbose && !is.null(X_full)) {
    message(sprintf("Constructed X_full_design with %d timepoints and %d regressors.", 
                    nrow(X_full), ncol(X_full)))
  } else if (verbose && is.null(X_full)) {
    message("Final X_full_design is NULL.")
  }
  X_full
}

# -------------------------------------------------------------------------
# Internal helper functions ------------------------------------------------
# -------------------------------------------------------------------------

#' @keywords internal
.ndx_validate_design_inputs <- function(run_idx,
                                        motion_params,
                                        rpca_components,
                                        spectral_sines,
                                        events = NULL) {

  if (length(run_idx) == 0) {
    stop("run_idx implies one or more runs have zero or negative length.")
  }

  if (any(is.na(run_idx) | !is.finite(run_idx))) {
    stop("run_idx contains NA or non-finite values.")
  }

  if (any(run_idx != floor(run_idx))) {
    stop("run_idx must contain integer values only.")
  }

  unique_vals   <- sort(unique(run_idx))
  run_mapping   <- setNames(seq_along(unique_vals), as.character(unique_vals))
  run_idx_mapped <- as.integer(run_mapping[as.character(run_idx)])
  unique_runs   <- sort(unique(run_idx_mapped))
  run_lengths   <- as.numeric(table(factor(run_idx_mapped, levels = unique_runs)))

  if (length(run_lengths) == 0 || sum(run_lengths) == 0) {
    stop("run_idx implies one or more runs have zero or negative length.")
  }

  if (any(run_lengths <= 0)) {
    stop("run_idx implies one or more runs have zero or negative length.")
  }
  total_tp <- sum(run_lengths)

  check_rows <- function(mat, name) {
    if (!is.null(mat) && nrow(mat) != total_tp) {
      stop(sprintf("Row mismatch: %s has %d rows, expected %d based on run_idx.",
                   name, nrow(mat), total_tp))
    }
  }
  check_rows(motion_params,   "motion_params")
  check_rows(rpca_components, "rpca_components")
  check_rows(spectral_sines,  "spectral_sines")

  if (!is.null(events)) {
    if (!is.data.frame(events) || !("blockids" %in% names(events))) {
      stop("events must be a data.frame containing a 'blockids' column for validation.")
    }
    invalid_blocks <- setdiff(unique(events$blockids), unique_runs)
    if (length(invalid_blocks) > 0) {
      stop(sprintf("events$blockids contains values not present in run_idx: %s",
                   paste(invalid_blocks, collapse = ", ")))
    }
  }

  list(unique_runs    = unique_runs,
       run_lengths    = run_lengths,
       total_tp       = total_tp,
       run_idx_mapped = run_idx_mapped,
       run_mapping    = run_mapping)

}

#' @keywords internal
.ndx_generate_task_regressors <- function(estimated_hrfs,
                                          events,
                                          sf, # sampling_frame from fmrireg
                                          TR,
                                          verbose = TRUE) {
  if (verbose) message("  Generating task regressors...") # Kept one high-level message if verbose
  
  if (is.null(estimated_hrfs) || !tibble::is_tibble(estimated_hrfs) || nrow(estimated_hrfs) == 0) {
    if (verbose) message("    No estimated_hrfs provided; skipping task regressor generation.")
    return(NULL)
  }
  if (!all(c("condition", "hrf_estimate") %in% names(estimated_hrfs))) {
    stop("estimated_hrfs tibble must contain 'condition' and 'hrf_estimate' columns.")
  }

  if (!is.null(events) && nrow(events) > 0) {
    if (!("blockids" %in% names(events))) {
      stop("events data frame must contain a 'blockids' column.")
    }
    n_runs <- length(sf$blocklens)
    if (any(is.na(events$blockids) | !is.finite(events$blockids) |
            events$blockids != floor(events$blockids))) {
      stop("events$blockids must be finite integers.")
    }
    if (any(!(events$blockids %in% seq_len(n_runs)))) {
      stop("events$blockids values outside range of run indices.")
    }
  }

  task_designs_list <- list()
  total_timepoints <- sum(sf$blocklens)

  for (i in 1:nrow(estimated_hrfs)) {
    cond_name  <- estimated_hrfs$condition[i]
    hrf_coeffs <- estimated_hrfs$hrf_estimate[[i]]
    
    if (is.null(hrf_coeffs) || !is.numeric(hrf_coeffs) || length(hrf_coeffs) == 0) {
      if (verbose) message(sprintf("    Skipping task regressor for condition '%s' due to NULL/empty or non-numeric HRF coefficients.", cond_name))
      next
    }
    num_fir <- length(hrf_coeffs)

    current_condition_events <- events[events$condition == cond_name, , drop=FALSE]
    if (nrow(current_condition_events) == 0) {
        if(verbose) message(sprintf("    No events found for condition '%s'. Skipping task regressor.", cond_name))
        next
    }

    X_fir_basis_cond <- matrix(0, nrow = total_timepoints, ncol = num_fir)
    fir_colnames <- paste0(make.names(cond_name), "_FIR", 0:(num_fir-1))
    colnames(X_fir_basis_cond) <- fir_colnames

    current_run_offset <- 0
    for (run_idx_val in seq_along(sf$blocklens)) {
      run_length_tps <- sf$blocklens[run_idx_val]
      run_start_tp_global <- current_run_offset + 1
      run_end_tp_global <- current_run_offset + run_length_tps
      run_events <- current_condition_events[current_condition_events$blockids == run_idx_val, ]

      if (nrow(run_events) > 0) {
        onset_tp_global <- floor(run_events$onsets / TR) + run_start_tp_global
        idx_matrix <- outer(onset_tp_global, 0:(num_fir - 1), "+")
        valid_mask <- idx_matrix >= run_start_tp_global &
                      idx_matrix <= run_end_tp_global &
                      idx_matrix <= total_timepoints
        if (any(valid_mask)) {
          row_idx <- idx_matrix[valid_mask]
          col_idx <- matrix(rep(seq_len(num_fir), each = length(onset_tp_global)),
                            nrow = length(onset_tp_global))[valid_mask]
          X_fir_basis_cond[cbind(row_idx, col_idx)] <- 1
        }
      }
      current_run_offset <- current_run_offset + run_length_tps
    }
    
    if (ncol(X_fir_basis_cond) == length(hrf_coeffs)) {
      X_task_cond <- X_fir_basis_cond %*% hrf_coeffs 
      colname_task <- paste0("task_", make.names(cond_name))
      colnames(X_task_cond) <- colname_task
      task_designs_list[[colname_task]] <- X_task_cond
    } else {
      if (verbose) message(sprintf("    Warning: Manual FIR basis for '%s' had %d columns, but hrf_coeffs length is %d. Skipping task regressor.", 
                                   cond_name, ncol(X_fir_basis_cond), length(hrf_coeffs)))
    }
  }
  
  if (length(task_designs_list) > 0) {
    final_task_matrix <- do.call(cbind, task_designs_list)
    return(final_task_matrix)
  } else {
    if (verbose) message("    No task regressors generated after iterating through conditions.")
    return(NULL)
  }
}

#' @keywords internal
.ndx_generate_nuisance_regressors <- function(motion_params,
                                              rpca_components,
                                              spectral_sines,
                                              verbose = TRUE) {
  if (verbose) message("  Generating nuisance regressors...") # Kept one high-level message
  nuis <- list()

  if ((is.null(motion_params) || ncol(motion_params) == 0) &&
      (is.null(rpca_components) || ncol(rpca_components) == 0) &&
      (is.null(spectral_sines) || ncol(spectral_sines) == 0)) {
    if (verbose) message("    No nuisance regressors provided.")
    return(NULL)
  }

  if (!is.null(motion_params) && ncol(motion_params) > 0) {
    nuis$motion <- as.matrix(motion_params)
    if (is.null(colnames(nuis$motion))) {
      colnames(nuis$motion) <- paste0("motion_param_", seq_len(ncol(nuis$motion)))
    }
  } else if (verbose && !is.null(motion_params)) {
    message("    motion_params provided but has 0 columns or is NULL.")
  }
  
  if (!is.null(rpca_components) && ncol(rpca_components) > 0) {
    nuis$rpca <- as.matrix(rpca_components)
    colnames(nuis$rpca) <- paste0("rpca_comp_", seq_len(ncol(nuis$rpca)))
  } else if (verbose && !is.null(rpca_components)) {
    message("    rpca_components provided but has 0 columns or is NULL.")
  }

  if (!is.null(spectral_sines) && ncol(spectral_sines) > 0) {
    nuis$spectral <- as.matrix(spectral_sines)
    if (is.null(colnames(nuis$spectral))) { 
        colnames(nuis$spectral) <- paste0("spectral_comp_", seq_len(ncol(nuis$spectral)))
    }
  } else if (verbose && !is.null(spectral_sines)) {
      message("    spectral_sines provided but has 0 columns or is NULL.")
  }
  nuis
}

#' @keywords internal
.ndx_generate_baseline_regressors <- function(run_idx,
                                             sf,
                                             poly_degree,
                                             verbose = TRUE) {
  if (verbose) message(sprintf("  Generating baseline regressors (poly_degree: %s)...", ifelse(is.null(poly_degree), "NULL", poly_degree))) # Kept one high-level
  if (any(is.na(run_idx) | !is.finite(run_idx) | run_idx != floor(run_idx))) {
    stop("run_idx for baseline must be finite integer values.")
  }
  unique_runs <- sort(unique(run_idx))
  total_tp    <- length(run_idx)
  baseline    <- list()
  has_poly0   <- FALSE

  if (!is.null(poly_degree) && is.numeric(poly_degree) && !is.na(poly_degree) && poly_degree >= 0) {
    if (poly_degree == 0) {
      if (verbose) message(sprintf("    Adding Legendre polynomial degree %d (intercept only).", poly_degree))
      bm_poly <- fmrireg::baseline_model(basis = "constant", intercept = "global", sframe = sf)
      X_poly  <- fmrireg::design_matrix(bm_poly)
      if (!is.null(X_poly) && ncol(X_poly) == 1) { 
        colnames(X_poly) <- "poly0"
        baseline$poly    <- X_poly
        has_poly0        <- TRUE
      } 
    } else { 
      if (verbose) message(sprintf("    Adding Legendre polynomials up to degree %d per run.", poly_degree))
      bm_poly <- fmrireg::baseline_model(basis = "poly", degree = poly_degree,
                                         intercept = "global", sframe = sf)
      X_poly  <- fmrireg::design_matrix(bm_poly)
      if (!is.null(X_poly) && ncol(X_poly) > 0) {
        colnames(X_poly) <- paste0("poly", 0:(ncol(X_poly)-1))
        baseline$poly    <- X_poly
        if ("poly0" %in% colnames(X_poly)) { 
            has_poly0 <- TRUE
        }
      } 
    }
  }

  if (length(unique_runs) > 1) {
    bm_run <- fmrireg::baseline_model(basis = "constant", degree = 0,
                                      intercept = "runwise", sframe = sf)
    X_run  <- fmrireg::design_matrix(bm_run)
    if (!is.null(X_run) && ncol(X_run) > 0) {
      if (has_poly0) {
        X_run <- X_run[, -1, drop = FALSE]
        unique_runs <- unique_runs[-1]
      }
      if (ncol(X_run) > 0) {
        colnames(X_run) <- paste0("run_intercept_", unique_runs)
        baseline$run_intercepts <- X_run
      }
    }
  } else if (length(unique_runs) == 1 && !has_poly0) {
    baseline$intercept <- matrix(1, nrow = total_tp, ncol = 1,
                                 dimnames = list(NULL, "intercept"))
    if (verbose) message("    Added global intercept for single run.")
  }

  if (length(baseline) > 0) {
    final_baseline_matrix <- as.matrix(do.call(cbind, unname(baseline)))
    return(final_baseline_matrix)
  } else {
    if (verbose) message("    No baseline regressors generated.")
    return(NULL)
  }
}

#' @keywords internal
#' @param zero_var_epsilon Numeric tolerance for detecting near-zero variance columns.
.ndx_combine_regressors <- function(reg_list,
                                   total_tp,
                                   verbose = TRUE,
                                   drop_zero_variance = FALSE,
                                   zero_var_epsilon = 1e-8) {
  valid_components <- list()
  for (name in names(reg_list)) {
    comp <- reg_list[[name]]
    is_comp_valid <- !is.null(comp) && is.matrix(comp) && ncol(comp) > 0
    if (is_comp_valid) {
      if (nrow(comp) != total_tp) {
        stop(sprintf("Component '%s' has %d rows, expected %d.", name, nrow(comp), total_tp))
      }
      valid_components[[name]] <- comp 
    } 
  }
  
  if (length(valid_components) == 0) {
    if (verbose) message("  No valid regressor components found to combine. Returning NULL.")
    return(NULL)
  }

  X_combined <- do.call(cbind, valid_components) 

  if (is.null(X_combined)) { 
      if (verbose) message("  Result of do.call(cbind,...) was NULL. Returning NULL.")
      return(NULL)
  }
  
  final_colnames <- unlist(lapply(valid_components, colnames))
  
  if (length(final_colnames) == ncol(X_combined)) {
      colnames(X_combined) <- make.names(final_colnames, unique = TRUE) 
  } else {
      warning("Mismatch between expected number of colnames and actual columns after cbind. Using default names.")
      colnames(X_combined) <- make.names(paste0("col", seq_len(ncol(X_combined))), unique = TRUE)
  }
  
   if (nrow(X_combined) != total_tp) {
     warning(sprintf("Combined design matrix has %d rows, expected %d. Critical error. Setting to NULL.",
                     nrow(X_combined), total_tp))
     return(NULL)
   }

  if (!is.null(X_combined) && ncol(X_combined) > 0) {
    col_vars <- apply(X_combined, 2, stats::var, na.rm = TRUE)
    potential_zero_var_cols <- colnames(X_combined)[col_vars < zero_var_epsilon]
    intercept_patterns <- "^(poly0|intercept|run_intercept_)"
    zero_var_cols_to_warn <- potential_zero_var_cols[!grepl(intercept_patterns, potential_zero_var_cols, ignore.case = TRUE)]
    if (length(zero_var_cols_to_warn) > 0 && verbose) {
      message(sprintf("  Warning: The following non-intercept columns in the final design matrix have near-zero variance: %s",
                      paste(zero_var_cols_to_warn, collapse=", ")))
    }
    if (drop_zero_variance && length(zero_var_cols_to_warn) > 0) {
      X_combined <- X_combined[, !(colnames(X_combined) %in% zero_var_cols_to_warn), drop = FALSE]
      if (verbose) message(sprintf("  Dropped %d near-zero variance column(s).", length(zero_var_cols_to_warn)))
    }
  }
  return(X_combined)
}
