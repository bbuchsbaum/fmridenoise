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
#' @param TR Repetition time, in seconds.
#' @param poly_degree Degree of Legendre polynomial baseline (per run).
#' @param verbose Logical flag for verbose console output.
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
                                    verbose = TRUE) {

  if (verbose) message("[ndx_build_design_matrix] Top level entry.")
  # 1. Basic validation ----------------------------------------------------
  info <- .ndx_validate_design_inputs(run_idx, motion_params, rpca_components, spectral_sines)
  if (verbose) message(sprintf("[ndx_build_design_matrix] Validation done. Total TP: %d", info$total_tp))
  sf   <- fmrireg::sampling_frame(blocklens = info$run_lengths, TR = TR)

  # 2. Task regressors -----------------------------------------------------
  if (verbose) message("[ndx_build_design_matrix] Calling .ndx_generate_task_regressors...")
  task_mat <- .ndx_generate_task_regressors(estimated_hrfs, events, sf, TR, verbose)
  if (verbose) {
    message("[ndx_build_design_matrix] Received task_mat:")
    if (is.null(task_mat)) message("  task_mat is NULL")
    else message(sprintf("  dim(task_mat): %s, colnames: %s", paste(dim(task_mat), collapse="x"), paste(colnames(task_mat), collapse=", ")))
  }

  # 3. Nuisance regressors -------------------------------------------------
  if (verbose) message("[ndx_build_design_matrix] Calling .ndx_generate_nuisance_regressors...")
  nuisance_list_components <- .ndx_generate_nuisance_regressors(motion_params, rpca_components, spectral_sines, verbose)
  if (verbose) {
    message("[ndx_build_design_matrix] Received nuisance_list_components:")
    if (length(nuisance_list_components) == 0) message("  nuisance_list_components is empty.")
    else for (n_name in names(nuisance_list_components)) {
      if (is.null(nuisance_list_components[[n_name]])) message(sprintf("  nuisance_list_components$%s is NULL", n_name))
      else message(sprintf("  nuisance_list_components$%s: dim: %s, colnames: %s", n_name, paste(dim(nuisance_list_components[[n_name]]), collapse="x"), paste(colnames(nuisance_list_components[[n_name]]), collapse=", ")))
    }
  }
  
  # 4. Baseline regressors -------------------------------------------------
  if (verbose) message("[ndx_build_design_matrix] Calling .ndx_generate_baseline_regressors...")
  baseline_mat <- .ndx_generate_baseline_regressors(run_idx, sf, poly_degree, verbose)
   if (verbose) {
    message("[ndx_build_design_matrix] Received baseline_mat:")
    if (is.null(baseline_mat)) message("  baseline_mat is NULL")
    else message(sprintf("  dim(baseline_mat): %s, colnames: %s", paste(dim(baseline_mat), collapse="x"), paste(colnames(baseline_mat), collapse=", ")))
  }

  # 5. Combine -------------------------------------------------------------
  regressor_list <- c(list(task = task_mat), nuisance_list_components, list(baseline = baseline_mat))
  if (verbose) {
    message("[ndx_build_design_matrix] Assembled regressor_list for combination:")
    if (length(regressor_list) == 0) message("  regressor_list is empty.")
    else for (r_name in names(regressor_list)) {
      if (is.null(regressor_list[[r_name]])) message(sprintf("  regressor_list$%s is NULL", r_name))
      else message(sprintf("  regressor_list$%s: dim: %s, colnames: %s", r_name, paste(dim(regressor_list[[r_name]]), collapse="x"), paste(colnames(regressor_list[[r_name]]), collapse=", ")))
    }
  }
  
  if (verbose) message("[ndx_build_design_matrix] Calling .ndx_combine_regressors...")
  X_full <- .ndx_combine_regressors(regressor_list, info$total_tp, verbose)

  if (verbose && !is.null(X_full)) {
    message(sprintf("[ndx_build_design_matrix] Final X_full: dim: %s, colnames: %s", 
                    paste(dim(X_full), collapse="x"), paste(colnames(X_full), collapse=", ")))
  } else if (verbose && is.null(X_full)) {
    message("[ndx_build_design_matrix] Final X_full is NULL.")
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
                                        spectral_sines) {
  # message("[.ndx_validate_design_inputs] Entry.") # Optional: too noisy if called often
  unique_runs <- sort(unique(run_idx))
  run_lengths <- as.numeric(table(factor(run_idx, levels = unique_runs)))

  # Added check for empty/zero-sum run_lengths to match test expectation
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

  list(unique_runs   = unique_runs,
       run_lengths   = run_lengths,
       total_tp      = total_tp)
}

#' @keywords internal
.ndx_generate_task_regressors <- function(estimated_hrfs,
                                          events,
                                          sf, # sampling_frame from fmrireg
                                          TR,
                                          verbose = TRUE) {
  if (verbose) message("[.ndx_generate_task_regressors] Entry (Manual FIR basis construction).")
  # Explicitly load fmrireg to be absolutely sure its functions are available in all contexts
  # This is a long shot, as devtools::load_all() should handle this.
  if (!("fmrireg" %in% .packages())) { # Only load if not already loaded by a direct library() call somewhere
      if (requireNamespace("fmrireg", quietly = TRUE)) {
          suppressPackageStartupMessages(library("fmrireg", character.only = TRUE)) # Corrected: package name as string
          if (verbose) message("  Dynamically loaded fmrireg namespace in .ndx_generate_task_regressors")
      } else {
          stop("fmrireg package is required but not found.")
      }
  }
  
  if (is.null(estimated_hrfs) || !tibble::is_tibble(estimated_hrfs) || nrow(estimated_hrfs) == 0) {
    if (verbose) message("  No estimated_hrfs provided; skipping task regressor generation.")
    return(NULL)
  }
  if (!all(c("condition", "hrf_estimate") %in% names(estimated_hrfs))) {
    stop("estimated_hrfs tibble must contain 'condition' and 'hrf_estimate' columns.")
  }

  task_designs_list <- list()
  total_timepoints <- sum(sf$blocklens)

  if (verbose) message(sprintf("  Generating task regressors for %d conditions. Total TPs: %d, TR: %.2f", nrow(estimated_hrfs), total_timepoints, TR))

  for (i in 1:nrow(estimated_hrfs)) {
    cond_name  <- estimated_hrfs$condition[i]
    hrf_coeffs <- estimated_hrfs$hrf_estimate[[i]]
    
    if (is.null(hrf_coeffs) || !is.numeric(hrf_coeffs) || length(hrf_coeffs) == 0) {
      if (verbose) message(sprintf("    Skipping task regressor for condition '%s' due to NULL/empty or non-numeric HRF coefficients.", cond_name))
      next
    }
    num_fir <- length(hrf_coeffs)
    if (verbose) message(sprintf("    Condition: %s, num_fir_taps: %d", cond_name, num_fir))

    current_condition_events <- events[events$condition == cond_name, , drop=FALSE]
    if (nrow(current_condition_events) == 0) {
        if(verbose) message(sprintf("    No events found for condition '%s'. Skipping task regressor.", cond_name))
        next
    }

    # Manual FIR basis matrix construction for this condition
    # Matrix will be total_timepoints x num_fir
    X_fir_basis_cond <- matrix(0, nrow = total_timepoints, ncol = num_fir)
    # Colnames for the basis (e.g., TaskA_FIR0, TaskA_FIR1, ...)
    fir_colnames <- paste0(make.names(cond_name), "_FIR", 0:(num_fir-1))
    colnames(X_fir_basis_cond) <- fir_colnames

    current_run_offset <- 0 # To handle timepoints across concatenated runs
    for (run_idx_val in 1:length(sf$blocklens)) {
      run_length_tps <- sf$blocklens[run_idx_val]
      run_start_tp_global <- current_run_offset + 1
      run_end_tp_global <- current_run_offset + run_length_tps
      
      # Get events for this condition that fall within this run
      # Note: event onsets are in seconds. sf$TR is seconds.
      # blockids in events_dm are 1-based and should align with run_idx_val
      run_events <- current_condition_events[current_condition_events$blockids == run_idx_val, ]

      if (nrow(run_events) > 0) {
        for (ev_idx in 1:nrow(run_events)) {
          event_onset_sec <- run_events$onsets[ev_idx]
          # Convert event onset from seconds to 0-indexed TRs *relative to the start of the current run*
          # Then add run_start_tp_global (1-indexed) and subtract 1 to get global 0-indexed TR
          event_onset_tr_global_0idx <- floor(event_onset_sec / TR) + (run_start_tp_global - 1)

          for (tap_idx in 0:(num_fir - 1)) { # 0-indexed tap
            target_timepoint_0idx <- event_onset_tr_global_0idx + tap_idx
            target_timepoint_1idx <- target_timepoint_0idx + 1

            if (target_timepoint_1idx >= run_start_tp_global && target_timepoint_1idx <= run_end_tp_global && target_timepoint_1idx <= total_timepoints) {
              # Ensure the target timepoint is within the bounds of the current run and overall duration
              X_fir_basis_cond[target_timepoint_1idx, tap_idx + 1] <- 1
            }
          }
        }
      }
      current_run_offset <- current_run_offset + run_length_tps
    }
    
    # Convolve FIR basis with HRF coefficients
    # hrf_coeffs should be a vector of length num_fir
    if (ncol(X_fir_basis_cond) == length(hrf_coeffs)) {
      X_task_cond <- X_fir_basis_cond %*% hrf_coeffs 
      # X_task_cond will be a single column matrix (total_timepoints x 1)
      colname_task <- paste0("task_", make.names(cond_name))
      colnames(X_task_cond) <- colname_task
      task_designs_list[[colname_task]] <- X_task_cond
      if (verbose) message(sprintf("      Generated task regressor '%s' using manual FIR. Dim: %s", colname_task, paste(dim(X_task_cond), collapse="x")))
    } else {
      if (verbose) message(sprintf("    Warning: Manual FIR basis for '%s' had %d columns, but hrf_coeffs length is %d. Skipping task regressor.", 
                                   cond_name, ncol(X_fir_basis_cond), length(hrf_coeffs)))
    }
  }

  if (verbose) message(sprintf("  [.ndx_generate_task_regressors] task_designs_list has %d elements before cbind.", length(task_designs_list)))
  
  if (length(task_designs_list) > 0) {
    final_task_matrix <- do.call(cbind, task_designs_list)
    if (verbose) {
        message(sprintf("  [.ndx_generate_task_regressors] Returning task matrix. dim: %s, colnames: %s", 
                        paste(dim(final_task_matrix), collapse="x"), paste(colnames(final_task_matrix), collapse=", ")))
    }
    return(final_task_matrix)
  } else {
    if (verbose) message("  [.ndx_generate_task_regressors] Returning NULL as no task regressors generated.")
    return(NULL)
  }
}

#' @keywords internal
.ndx_generate_nuisance_regressors <- function(motion_params,
                                              rpca_components,
                                              spectral_sines,
                                              verbose = TRUE) {
  if (verbose) message("[.ndx_generate_nuisance_regressors] Entry.")
  nuis <- list()

  if (!is.null(motion_params) && ncol(motion_params) > 0) {
    nuis$motion <- as.matrix(motion_params)
    if (is.null(colnames(nuis$motion))) {
      colnames(nuis$motion) <- paste0("motion_param_", seq_len(ncol(nuis$motion)))
    }
    if (verbose) message(sprintf("  Added motion_params. dim: %s, colnames: %s", paste(dim(nuis$motion), collapse="x"), paste(colnames(nuis$motion), collapse=", ")))
  } else if (verbose && !is.null(motion_params)) {
    message("  motion_params provided but has 0 columns or is NULL.")
  }
  
  if (!is.null(rpca_components) && ncol(rpca_components) > 0) {
    nuis$rpca <- as.matrix(rpca_components)
    colnames(nuis$rpca) <- paste0("rpca_comp_", seq_len(ncol(nuis$rpca)))
    if (verbose) message(sprintf("  Added rpca_components. dim: %s, colnames: %s", paste(dim(nuis$rpca), collapse="x"), paste(colnames(nuis$rpca), collapse=", ")))
  } else if (verbose && !is.null(rpca_components)) {
    message("  rpca_components provided but has 0 columns or is NULL.")
  }

  if (!is.null(spectral_sines) && ncol(spectral_sines) > 0) {
    nuis$spectral <- as.matrix(spectral_sines)
    if (is.null(colnames(nuis$spectral))) { 
        colnames(nuis$spectral) <- paste0("spectral_comp_", seq_len(ncol(nuis$spectral)))
    }
    if (verbose) message(sprintf("  Added spectral_sines. dim: %s, colnames: %s", paste(dim(nuis$spectral), collapse="x"), paste(colnames(nuis$spectral), collapse=", ")))
  } else if (verbose && !is.null(spectral_sines)) {
      message("  spectral_sines provided but has 0 columns or is NULL.")
  }
  
  if (verbose) message(sprintf("[.ndx_generate_nuisance_regressors] Returning nuisance list with %d elements. Names: %s", length(nuis), paste(names(nuis),collapse=", ")))
  nuis
}

#' @keywords internal
.ndx_generate_baseline_regressors <- function(run_idx,
                                             sf,
                                             poly_degree,
                                             verbose = TRUE) {
  if (verbose) message(sprintf("[.ndx_generate_baseline_regressors] Entry. poly_degree: %s", ifelse(is.null(poly_degree), "NULL", poly_degree)))
  unique_runs <- sort(unique(run_idx))
  total_tp    <- length(run_idx)
  baseline    <- list()
  has_poly0   <- FALSE

  # Legendre polynomials ---------------------------------------------------
  if (!is.null(poly_degree) && is.numeric(poly_degree) && !is.na(poly_degree) && poly_degree >= 0) {
    if (poly_degree == 0) {
      # Handle intercept-only case (poly0) using basis = "constant"
      if (verbose) message(sprintf("  Adding Legendre polynomial degree %d (intercept only).", poly_degree))
      bm_poly <- fmrireg::baseline_model(basis = "constant", intercept = "global", sframe = sf)
      X_poly  <- fmrireg::design_matrix(bm_poly)
      if (!is.null(X_poly) && ncol(X_poly) == 1) { # Expect 1 column for intercept
        colnames(X_poly) <- "poly0"
        baseline$poly    <- X_poly
        has_poly0        <- TRUE
        if (verbose) message(sprintf("    Added 1 polynomial baseline regressor (poly0)."))
      } else if (verbose) {
          message(sprintf("    Warning: Expected 1 column for poly0 (intercept), got %s. Not adding poly0.", if(is.null(X_poly)) "NULL" else as.character(ncol(X_poly))))
      }
    } else { # poly_degree > 0
      if (verbose) message(sprintf("  Adding Legendre polynomials up to degree %d per run.", poly_degree))
      bm_poly <- fmrireg::baseline_model(basis = "poly", degree = poly_degree,
                                         intercept = "global", sframe = sf)
      X_poly  <- fmrireg::design_matrix(bm_poly)
      if (!is.null(X_poly) && ncol(X_poly) > 0) {
        # Name them poly0, poly1, ... up to ncol-1
        colnames(X_poly) <- paste0("poly", 0:(ncol(X_poly)-1))
        baseline$poly    <- X_poly
        if ("poly0" %in% colnames(X_poly)) { # Check if an intercept term is present
            has_poly0 <- TRUE
        }
        if (verbose) message(sprintf("    Added %d polynomial baseline regressors from fmrireg. Colnames: %s", ncol(X_poly), paste(colnames(X_poly), collapse=", ")))
      } else if (verbose) {
          message(sprintf("    Polynomial basis generation (degree %d) resulted in NULL or 0 columns.", poly_degree))
      }
    }
  }

  # Run intercepts ---------------------------------------------------------
  if (length(unique_runs) > 1) {
    bm_run <- fmrireg::baseline_model(basis = "constant", degree = 0,
                                      intercept = "runwise", sframe = sf)
    X_run  <- fmrireg::design_matrix(bm_run)
    if (!is.null(X_run) && ncol(X_run) > 0) {
      if (has_poly0) {
        # Drop the first column to avoid redundant intercept
        X_run <- X_run[, -1, drop = FALSE]
        unique_runs <- unique_runs[-1]
      }
      if (ncol(X_run) > 0) {
        colnames(X_run) <- paste0("run_intercept_", unique_runs)
        baseline$run_intercepts <- X_run
        if (verbose) message(sprintf("  Added %d run intercept regressors.", ncol(X_run)))
      }
    }
  } else if (length(unique_runs) == 1 && !has_poly0) {
    baseline$intercept <- matrix(1, nrow = total_tp, ncol = 1,
                                 dimnames = list(NULL, "intercept"))
    if (verbose) message("  Added global intercept.")
  }

  if (verbose) message(sprintf("  [.ndx_generate_baseline_regressors] baseline list has %d elements before cbind. Names: %s", length(baseline), paste(names(baseline),collapse=", ")))

  if (length(baseline) > 0) {
    final_baseline_matrix <- as.matrix(do.call(cbind, unname(baseline))) # Use unname() to avoid list name prefixes
    if (verbose) {
        message(sprintf("  [.ndx_generate_baseline_regressors] Returning baseline matrix. dim: %s, colnames: %s", 
                        paste(dim(final_baseline_matrix), collapse="x"), paste(colnames(final_baseline_matrix), collapse=", ")))
    }
    return(final_baseline_matrix)
  } else {
    if (verbose) message("  [.ndx_generate_baseline_regressors] Returning NULL as no baseline regressors generated.")
    return(NULL)
  }
}

#' @keywords internal
.ndx_combine_regressors <- function(reg_list,
                                   total_tp,
                                   verbose = TRUE) {
  if (verbose) {
      message("[.ndx_combine_regressors] Entry. Input reg_list:")
      if (length(reg_list) == 0) message("  Input reg_list is empty.")
      else for (r_name in names(reg_list)) {
        if (is.null(reg_list[[r_name]])) message(sprintf("  Input reg_list$%s is NULL", r_name))
        else message(sprintf("  Input reg_list$%s: dim: %s, colnames: %s", r_name, paste(dim(reg_list[[r_name]]), collapse="x"), paste(colnames(reg_list[[r_name]]), collapse=", ")))
      }
  }
  valid_components <- list()

  for (name in names(reg_list)) {
    comp <- reg_list[[name]]
    is_comp_valid <- !is.null(comp) && is.matrix(comp) && ncol(comp) > 0
    
    if (is_comp_valid) {
      if (nrow(comp) != total_tp) {
        stop(sprintf("Component '%s' has %d rows, expected %d.", name, nrow(comp), total_tp))
      }
      valid_components[[name]] <- comp 
      if (verbose) message(sprintf("  Added valid component: %s. dim: %s", name, paste(dim(comp),collapse="x")))
    } else { 
        if (verbose) {
            reason <- c()
            if (is.null(comp)) reason <- c(reason, "is NULL")
            # Check !is.matrix only if not NULL, otherwise is.matrix(NULL) is FALSE and might be confusing.
            if (!is.null(comp) && !is.matrix(comp)) reason <- c(reason, "not a matrix")
            # Check ncol only if it is a matrix.
            if (!is.null(comp) && is.matrix(comp) && ncol(comp) == 0) reason <- c(reason, "0 cols")
            # Fallback if no specific reason caught but comp is not valid
            if (length(reason) == 0 && !is_comp_valid) reason <- c(reason, "failed validation (unexpected)")
            message(sprintf("  Skipping component '%s'. Reasons: %s.", name, paste(reason, collapse=", ")))
        }
    }
  }

  if (verbose) {
    message("[.ndx_combine_regressors] valid_components collected:")
    if (length(valid_components) == 0) message("  valid_components is empty.")
    else for (v_name in names(valid_components)) {
        message(sprintf("  valid_components$%s: dim: %s, colnames: %s", v_name, paste(dim(valid_components[[v_name]]), collapse="x"), paste(colnames(valid_components[[v_name]]), collapse=", ")))
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
  if (verbose) message(sprintf("  X_combined after cbind: dim: %s, initial colnames: %s", paste(dim(X_combined),collapse="x"), paste(colnames(X_combined),collapse=", ")))
  
  final_colnames <- unlist(lapply(valid_components, colnames))
  if (verbose) message(sprintf("  final_colnames from valid_components: %s", paste(final_colnames,collapse=", ")))
  
  if (length(final_colnames) == ncol(X_combined)) {
      # If we have the correct number of names, and they come from components that should already be well-formed,
      # apply them directly first, then ensure uniqueness with make.names if truly needed (e.g. if cbind somehow lost uniqueness).
      # For now, let's assume components provide unique enough names that make.names primarily handles syntax.
      colnames(X_combined) <- make.names(final_colnames, unique = TRUE) # Keep make.names for syntax and final uniqueness pass
  } else {
      warning("Mismatch between expected number of colnames and actual columns after cbind. Using default names.")
      colnames(X_combined) <- make.names(paste0("col", seq_len(ncol(X_combined))), unique = TRUE)
  }
  if (verbose) message(sprintf("  X_combined after make.names: dim: %s, final colnames: %s", paste(dim(X_combined),collapse="x"), paste(colnames(X_combined),collapse=", ")))
  
   if (nrow(X_combined) != total_tp) {
     warning(sprintf("Combined design matrix has %d rows, expected %d. Critical error. Setting to NULL.",
                     nrow(X_combined), total_tp))
     return(NULL)
   }

  # Check for zero-variance columns (excluding intercept-like terms)
  if (!is.null(X_combined) && ncol(X_combined) > 0) {
    col_vars <- apply(X_combined, 2, stats::var, na.rm = TRUE)
    # Identify columns with variance very close to zero
    # Exclude known intercept-like columns from this specific warning, as they are often intentionally constant.
    # A more robust check might involve looking at the type of regressor if that info was available.
    potential_zero_var_cols <- colnames(X_combined)[col_vars < (.Machine$double.eps * 100)] # Use a slightly larger epsilon
    
    # Filter out known intercept column names (case-insensitive)
    intercept_patterns = "^(poly0|intercept|run_intercept_)"
    zero_var_cols_to_warn <- potential_zero_var_cols[!grepl(intercept_patterns, potential_zero_var_cols, ignore.case = TRUE)]
    
    if (length(zero_var_cols_to_warn) > 0 && verbose) {
      message(sprintf("  Warning: The following non-intercept columns in the final design matrix have near-zero variance: %s", 
                      paste(zero_var_cols_to_warn, collapse=", ")))
    }
  }

  if (verbose) message("[.ndx_combine_regressors] Returning X_combined.")
  return(X_combined)
} 