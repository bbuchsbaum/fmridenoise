#' Run RPCA for a Single Run
#'
#' Internal helper used by `ndx_rpca_temporal_components_multirun`.
#' Performs RPCA on a single run's residual matrix and extracts the
#' voxel-space components.
#'
#' @param Er Numeric matrix of residuals for one run (Time x Voxels).
#' @param run_name Character run label used in messages.
#' @param opts List of options controlling the RPCA step.
#' @return List with elements `V_r`, `S_r_t`, `spike_mask` and `glitch_ratio`.
#' @keywords internal
.ndx_rpca_single_run <- function(Er, run_name, opts) {
  spike_mask <- if (nrow(Er) > 0) rep(FALSE, nrow(Er)) else logical(0)

  if (nrow(Er) == 0 || ncol(Er) == 0) {
    warning(sprintf("Residuals for run %s are empty (dims: %s). Skipping RPCA for this run.",
                    run_name, paste(dim(Er), collapse = "x")))
    return(list(V_r = NULL, S_r_t = NULL, spike_mask = spike_mask, glitch_ratio = NA))
  }

  Er_t <- t(Er)
  lambda_r <- if (opts$rpca_lambda_auto) {
    1 / sqrt(max(dim(Er_t)))
  } else {
    if (is.null(opts$rpca_lambda_fixed))
      stop("rpca_lambda_fixed must be provided if rpca_lambda_auto is FALSE.")
    opts$rpca_lambda_fixed
  }

  rpca_call_args <- list(
    M = Er_t,
    lambda = lambda_r,
    term.delta = opts$rpca_term_delta,
    max.iter = opts$rpca_max_iter,
    trace = opts$rpca_trace
  )

  if (!is.null(opts$rpca_mu)) {
    rpca_call_args$mu <- opts$rpca_mu
    message(sprintf("  Run %s: Using user-specified mu = %f for rpca.", run_name, opts$rpca_mu))
  } else {
    temp_mu_for_logging <- if (sum(abs(Er_t)) > 1e-9) {
      prod(dim(Er_t)) / (4 * sum(abs(Er_t)))
    } else {
      1.0
    }
    message(sprintf(
      "  Run %s: rpca_mu is NULL, rpca::rpca will use its internal default mu (approx for logging: %.2e).",
      run_name, temp_mu_for_logging
    ))
  }

  k_this_run <- min(opts$k_per_run_target, min(dim(Er_t)))
  if (k_this_run <= 0) {
    warning(sprintf(
      "Cannot perform RPCA for run %s: k_this_run (%d) is not positive after adjustment. Skipping.",
      run_name, k_this_run
    ))
    return(list(V_r = NULL, S_r_t = NULL, spike_mask = spike_mask, glitch_ratio = NA))
  }

  message(sprintf(
    "  Processing run %s (data: %d voxels x %d TRs, target k: %d, lambda: %.2e)",
    run_name, nrow(Er_t), ncol(Er_t), k_this_run, lambda_r
  ))

  rpca_res_r <- NULL
  tryCatch({
    rpca_res_r <- do.call(rpca, rpca_call_args)
  }, error = function(e) {
    warning(sprintf("rpca::rpca failed for run %s: %s", run_name, e$message))
    rpca_res_r <<- NULL
  })

  if (is.null(rpca_res_r) || is.null(rpca_res_r$L)) {
    warning(sprintf("Robust PCA failed to produce L for run %s (or rpca_res_r is NULL).", run_name))
    return(list(V_r = NULL, S_r_t = NULL, spike_mask = spike_mask, glitch_ratio = NA))
  }

  L_r_t <- rpca_res_r$L
  if (is.null(L_r_t) || min(dim(L_r_t)) == 0 || sum(abs(L_r_t)) < 1e-9) {
    warning(sprintf("RPCA for run %s yielded an empty or zero L component.", run_name))
    return(list(V_r = NULL, S_r_t = rpca_res_r$S, spike_mask = spike_mask, glitch_ratio = NA))
  }

  k_eff_svd_L <- min(k_this_run, nrow(L_r_t), ncol(L_r_t))
  if (k_eff_svd_L < k_this_run) {
    message(sprintf(
      "  Run %s: Effective k for SVD of L_r_t (%d) is less than k_this_run (%d) due to matrix dimensions.",
      run_name, k_eff_svd_L, k_this_run
    ))
  }

  svd_L_r_t <- NULL
  tryCatch({
    svd_L_r_t <- svd(L_r_t, nu = k_eff_svd_L, nv = 0)
  }, error = function(e) {
    warning(sprintf("SVD on L_r_t for run %s failed: %s", run_name, e$message))
    svd_L_r_t <<- NULL
  })

  if (is.null(svd_L_r_t) || is.null(svd_L_r_t$u) || ncol(svd_L_r_t$u) == 0) {
    warning(sprintf("SVD of L_r_t for run %s yielded no U components for V_r.", run_name))
    return(list(V_r = NULL, S_r_t = rpca_res_r$S, spike_mask = spike_mask, glitch_ratio = NA))
  }

  V_r <- svd_L_r_t$u
  S_r_t <- rpca_res_r$S

  glitch_ratio <- NA
  if (!is.null(S_r_t) && !is.null(L_r_t)) {
    energy_S_r <- sum(S_r_t^2)
    energy_L_r <- sum(L_r_t^2)
    if (energy_L_r > 1e-9) glitch_ratio <- energy_S_r / energy_L_r
  }

  if (!is.null(S_r_t) && length(S_r_t) > 0 && sum(abs(S_r_t), na.rm = TRUE) > 1e-9) {
    s_t_star_run <- apply(abs(S_r_t), 2, stats::median, na.rm = TRUE)
    mad_s_t_star <- stats::mad(s_t_star_run, constant = 1, na.rm = TRUE)
    if (!is.finite(mad_s_t_star) || mad_s_t_star < 1e-9) {
      thresh_val <- stats::quantile(s_t_star_run,
                                    probs = opts$rpca_spike_percentile_thresh,
                                    na.rm = TRUE)
    } else {
      thresh_val <- median(s_t_star_run, na.rm = TRUE) +
        opts$rpca_spike_mad_thresh * mad_s_t_star
    }
    spike_tmp <- s_t_star_run > thresh_val
    if (length(spike_tmp) == nrow(Er)) {
      spike_mask <- as.logical(spike_tmp)
    } else {
      warning(sprintf(
        "Run %s: Generated spike_TR_mask_run length (%d) did not match TRs for run (%d). Using all FALSEs.",
        run_name, length(spike_tmp), nrow(Er)
      ))
      spike_mask <- rep(FALSE, nrow(Er))
    }
  }

  list(V_r = V_r, S_r_t = S_r_t, spike_mask = spike_mask, glitch_ratio = glitch_ratio)
}

#' Merge Voxel Components Across Runs
#'
#' Handles the choice between `concat_svd` and `iterative` strategies when
#' forming the global voxel-space basis.
#'
#' @param V_list List of per-run V_r matrices (may contain `NULL`).
#' @param opts Options list controlling the merge strategy.
#' @param k_target Desired number of global components.
#' @return List with `V_global` and `singular_values`.
#' @keywords internal
.ndx_merge_voxel_components <- function(V_list, opts, k_target) {
  V_list_valid <- Filter(Negate(is.null), V_list)
  if (length(V_list_valid) == 0) {
    warning("RPCA failed for all runs, or no components were extracted. Cannot proceed.")
    return(list(V_global = NULL, singular_values = NULL))
  }
  message(sprintf("Successfully extracted initial voxel components (V_r) from %d runs.",
                  length(V_list_valid)))

  V_global <- NULL
  V_global_singular_values <- NULL

  if (opts$rpca_merge_strategy == "iterative") {
    message("Using Iterative Grassmann Averaging for V_global...")
    V_global <- .grassmann_merge_iterative(V_list_valid, k_target)
    if (is.null(V_global)) {
      warning("Iterative Grassmann Averaging failed to produce V_global.")
      return(list(V_global = NULL, singular_values = NULL))
    }
  } else if (opts$rpca_merge_strategy == "concat_svd") {
    message("Using Concatenate & SVD method for V_global...")
    V_all_concat_list <- Filter(function(x) !is.null(x) && ncol(x) > 0, V_list_valid)
    if (length(V_all_concat_list) == 0) {
      warning("No valid V_r components to concatenate for SVD V_global step.")
      return(list(V_global = NULL, singular_values = NULL))
    }
    V_all_concat <- do.call(cbind, V_all_concat_list)
    if (is.null(V_all_concat) || ncol(V_all_concat) == 0) {
      warning("Concatenation of V_r components for SVD V_global resulted in an empty matrix.")
      return(list(V_global = NULL, singular_values = NULL))
    }
    message(sprintf(
      "Performing SVD on concatenated V_r matrix (dims: %s) to find V_global (%d components)",
      paste(dim(V_all_concat), collapse = "x"), k_target
    ))
    k_for_svd_V_all <- min(k_target, ncol(V_all_concat), nrow(V_all_concat))
    if (k_for_svd_V_all <= 0) {
      warning(sprintf(
        "Cannot perform SVD on concatenated V_r: k_for_svd_V_all (%d) is not positive.",
        k_for_svd_V_all
      ))
      return(list(V_global = NULL, singular_values = NULL))
    }
    svd_V_all <- NULL
    tryCatch({
      svd_V_all <- svd(V_all_concat, nu = k_for_svd_V_all, nv = 0)
    }, error = function(e) {
      warning(paste("SVD on concatenated V_r components for V_global failed:", e$message))
      svd_V_all <<- NULL
    })
    if (is.null(svd_V_all) || is.null(svd_V_all$u) || ncol(svd_V_all$u) == 0) {
      warning("SVD on concatenated V_r for V_global failed to produce U vectors.")
      return(list(V_global = NULL, singular_values = NULL))
    }
    V_global <- svd_V_all$u
    V_global_singular_values <- svd_V_all$d
    message(sprintf("  V_global (concat_svd) obtained with %d components.", ncol(V_global)))
  } else {
    stop(sprintf("Invalid rpca_merge_strategy: '%s'. Choose 'concat_svd' or 'iterative'.",
                 opts$rpca_merge_strategy))
  }

  if (is.null(V_global) || ncol(V_global) == 0) {
    warning("V_global is NULL or has zero components after merging. Cannot proceed.")
    return(list(V_global = NULL, singular_values = NULL))
  }

  list(V_global = V_global, singular_values = V_global_singular_values)
}

#' Form Temporal Components From Residuals
#'
#' Multiplies each run's residuals by the global voxel-space basis and
#' concatenates the results.
#'
#' @param residuals_list List of residual matrices (Time x Voxels) for each run.
#' @param V_global Voxel-space basis matrix returned from merging.
#' @return Matrix of concatenated temporal components.
#' @keywords internal
.ndx_form_temporal_components <- function(residuals_list, V_global) {
  message("Forming run-specific temporal components C_r...")
  C_r_list <- lapply(residuals_list, function(E) {
    if (is.null(E) || nrow(E) == 0) {
      matrix(0, nrow = nrow(E), ncol = ncol(V_global))
    } else {
      E %*% V_global
    }
  })
  do.call(rbind, C_r_list)
}

