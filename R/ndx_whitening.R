#' AR(2) Pre-whitening for fMRI Data
#'
#' Estimates AR(2) model coefficients voxel-wise from provided residuals and applies
#' the whitening transformation to the fMRI data and the design matrix.
#'
#' @param Y_data A numeric matrix (timepoints x voxels) of fMRI data to be whitened.
#' @param X_design_full A numeric matrix (timepoints x regressors) representing the
#'   full design matrix to be whitened.
#' @param Y_residuals_for_AR_fit A numeric matrix (timepoints x voxels) of residuals
#'   used to estimate the AR(2) coefficients. This should typically be from a model
#'   that accounts for task effects and other known structured noise components.
#' @param order Integer, the order of the AR model. Defaults to 2 for AR(2).
#' @param global_ar_on_design Logical, if TRUE (default), a global AR model (averaged from successful
#'   voxel-wise fits) is used to whiten `X_design_full`. If FALSE, `X_design_full` is not whitened,
#'   and users should handle its whitening if necessary for downstream voxel-wise modeling.
#' @param verbose Logical, if TRUE (default) progress messages are printed. Set to FALSE to silence
#'   console output.
#' @param weights Optional numeric matrix of weights (timepoints x voxels) used
#'   for weighted AR coefficient estimation. If NULL, unweighted estimation is
#'   performed.
#' @param max_ar_failures_prop Numeric between 0 and 1. If the proportion of
#'   voxels with failed or unstable AR fits exceeds this threshold, whitening is
#'   skipped and a warning is issued. Default: 0.3.
#' @param global_ar_stat Character string specifying how to summarize voxel-wise
#'   AR coefficients when computing the global filter for `X_design_full`.
#'   Options are "mean" (default), "median", or "trimmed_mean" (10% trimmed
#'   mean).
#'
#' @return A list containing:
#'   - `Y_whitened`: The AR(order)-whitened fMRI data (timepoints x voxels).
#'   - `X_whitened`: The AR(order)-whitened design matrix (timepoints x regressors). If `global_ar_on_design` is FALSE,
#'     this will be the original `X_design_full`.
#'   - `AR_coeffs_voxelwise`: A matrix (voxels x order) of the estimated voxel-wise AR coefficients.
#'     Stable coefficients are returned; unstable ones are set to zero. The first column corresponds to phi_1, the second to phi_2, etc.
#'   - `AR_coeffs_global`: A numeric vector of length `order` representing the global AR coefficients used for
#'     whitening `X_design_full` (if `global_ar_on_design` is TRUE). NULL otherwise.
#'   - `var_innovations_voxelwise`: A numeric vector (voxels) of the estimated voxel-wise innovation variances
#'     (i.e., variance of residuals after whitening). NA for voxels where AR fit failed.
#'   - `na_mask`: A logical vector of length `n_timepoints`. TRUE for the initial `order` rows that are NA in
#'     `Y_whitened` and `X_whitened` (if applicable), FALSE otherwise.
#'
#' @details
#' The AR(order) model for a time series `y_t` is: 
#' `y_t = phi_1 * y_{t-1} + ... + phi_order * y_{t-order} + e_t`, where `e_t` is white noise.
#' The whitened value `y*_t` (i.e., the estimate of `e_t`) is `y_t - phi_1 * y_{t-1} - ... - phi_order * y_{t-order}`.
#' This transformation is applied to `Y_data` using voxel-specific coefficients. Unstable AR coefficients are set to zero (no whitening for that voxel).
#' If `global_ar_on_design` is TRUE, `X_design_full` is whitened using averaged AR coefficients (from stable fits).
#' Whitening is performed with `stats::filter(method = "convolution", sides = 1)` so that the
#' output represents the innovations of the AR model.
#' The first `order` rows of the whitened matrices will contain `NA`s due to the filter initialization.
#' The `na_mask` in the output identifies these rows.
#'
#' @import stats
#' @export
ndx_ar2_whitening <- function(Y_data, X_design_full, Y_residuals_for_AR_fit,
                              order = 2L, global_ar_on_design = TRUE,
                              verbose = TRUE, weights = NULL,
                              max_ar_failures_prop = 0.3,
                              global_ar_stat = c("mean", "median", "trimmed_mean")) {

  if (!is.matrix(Y_data) || !is.numeric(Y_data)) {
    stop("Y_data must be a numeric matrix.")
  }
  if (!is.matrix(X_design_full) || !is.numeric(X_design_full)) {
    stop("X_design_full must be a numeric matrix.")
  }
  if (!is.matrix(Y_residuals_for_AR_fit) || !is.numeric(Y_residuals_for_AR_fit)) {
    stop("Y_residuals_for_AR_fit must be a numeric matrix.")
  }
  if (nrow(Y_data) != nrow(X_design_full) || nrow(Y_data) != nrow(Y_residuals_for_AR_fit)) {
    stop("Y_data, X_design_full, and Y_residuals_for_AR_fit must have the same number of rows (timepoints).")
  }
  if (ncol(Y_data) != ncol(Y_residuals_for_AR_fit)) {
    stop("Y_data and Y_residuals_for_AR_fit must have the same number of columns (voxels).")
  }
  if (!is.numeric(order) || length(order) != 1 || !is.finite(order) || order < 1 || order %% 1 != 0) {
    stop("order must be a single positive integer.")
  }
  order <- as.integer(order)
  global_ar_stat <- match.arg(global_ar_stat)

  n_timepoints <- nrow(Y_data)
  n_voxels <- ncol(Y_data)

  if (n_timepoints <= order) {
    stop(sprintf("Number of timepoints (%d) must be greater than AR order (%d).", n_timepoints, order))
  }

  if (!is.null(weights)) {
    if (!is.matrix(weights) || !is.numeric(weights) ||
        nrow(weights) != n_timepoints || ncol(weights) != n_voxels) {
      stop("weights must be a numeric matrix with dimensions matching Y_data")
    }
    if (any(!is.finite(weights)) || any(weights < 0)) {
      stop("weights must contain finite, non-negative numbers")
    }
  }

  if (verbose) message(sprintf("Estimating AR(%d) coefficients for %d voxels...", order, n_voxels))

  fit_env <- new.env(parent = emptyenv())
  fit_env$idx <- 0L

  ar_results <- apply(Y_residuals_for_AR_fit, 2, function(res_col) {
    fit_env$idx <- fit_env$idx + 1L
    w_vec <- if (is.null(weights)) NULL else weights[, fit_env$idx]
    .fit_ar_single_voxel(res_col, order = order, weights_vec = w_vec)
  })

  if (is.vector(ar_results)) ar_results <- matrix(ar_results, nrow = order + 1)

  AR_coeffs_voxelwise <- t(ar_results[seq_len(order), , drop = FALSE])
  var_innovations_voxelwise <- as.numeric(ar_results[order + 1, ])
  num_phi_zeroed <- sum(rowSums(AR_coeffs_voxelwise != 0) == 0)

  if (verbose) message("AR coefficient estimation complete.")

  # Report on initially NA coefficients (which are now 0 if stability check also failed or fit was null)
  # num_na_initial <- sum(is.na(var_innovations_voxelwise)) # More direct count of actual fit failures before stability
  # if (num_na_initial > 0) {
  #     message(sprintf("AR fitting initially failed for %d/%d voxels (e.g., due to low variance).", num_na_initial, n_voxels))
  # }
  
  if (num_phi_zeroed > 0 && verbose) {
      message(sprintf("%d/%d voxels had AR coefficients set to zero (due to fit failure or instability). These will not be whitened.", num_phi_zeroed, n_voxels))
  }
  
  failure_prop <- if (n_voxels > 0) num_phi_zeroed / n_voxels else 0
  if (failure_prop > max_ar_failures_prop) {
      warning(sprintf(
        "Proportion of failed AR fits (%.2f) exceeds max_ar_failures_prop %.2f. Skipping whitening.",
        failure_prop, max_ar_failures_prop))
      return(list(
        Y_whitened = Y_data,
        X_whitened = X_design_full,
        AR_coeffs_voxelwise = AR_coeffs_voxelwise,
        AR_coeffs_global = NULL,
        var_innovations_voxelwise = var_innovations_voxelwise,
        na_mask = rep(FALSE, n_timepoints)
      ))
  }

  if (verbose) message("Applying voxel-specific AR filter to Y_data...")
  Y_whitened <- .apply_ar_filter_voxelwise(Y_data, AR_coeffs_voxelwise, order)
  if (verbose) message("Y_data whitening complete.")

  AR_coeffs_global <- NULL
  if (global_ar_on_design) {
    if (verbose) message("Calculating global AR coefficients for X_design_full...")
    valid_coeffs_for_global_avg <- AR_coeffs_voxelwise[rowSums(AR_coeffs_voxelwise != 0) > 0 & !is.na(var_innovations_voxelwise), , drop = FALSE]

    if (nrow(valid_coeffs_for_global_avg) > 0) {
      AR_coeffs_global <- switch(global_ar_stat,
                                 mean = colMeans(valid_coeffs_for_global_avg, na.rm = TRUE),
                                 median = apply(valid_coeffs_for_global_avg, 2, stats::median, na.rm = TRUE),
                                 trimmed_mean = apply(valid_coeffs_for_global_avg, 2, mean, trim = 0.1, na.rm = TRUE))
      if(any(is.na(AR_coeffs_global))) {
          warning("Global AR coefficients for design matrix contained NAs after averaging. Design matrix will not be whitened.")
          X_whitened <- X_design_full
          AR_coeffs_global <- NULL
      } else {
          # Stability check on averaged coefficients
          roots_global <- tryCatch(polyroot(c(1, -AR_coeffs_global)), error = function(e) NULL)
          if(is.null(roots_global) || any(is.na(roots_global)) || any(Mod(roots_global) <= 1.00001)) {
              warning("Averaged AR coefficients for design matrix are unstable. Design matrix will not be whitened.")
              X_whitened <- X_design_full
              AR_coeffs_global <- NULL
          } else {
              if (verbose) message(sprintf("Applying global AR filter (coeffs: %s) to X_design_full...", paste(round(AR_coeffs_global,3), collapse=", ")))
              X_whitened <- .apply_ar_filter_to_matrix_cols(X_design_full, AR_coeffs_global, order)
              if (verbose) message("X_design_full whitening complete.")
          }
      }
    } else {
      warning("No valid voxel-wise AR coefficients available to compute global AR model for design matrix. Design matrix will not be whitened.")
      X_whitened <- X_design_full
    }
  } else {
    if (verbose) message("Skipping whitening of X_design_full as per global_ar_on_design = FALSE.")
    X_whitened <- X_design_full
  }
  
  na_mask <- seq_len(n_timepoints) <= order
  
  return(list(
    Y_whitened = Y_whitened,
    X_whitened = X_whitened,
    AR_coeffs_voxelwise = AR_coeffs_voxelwise,
    AR_coeffs_global = AR_coeffs_global,
    var_innovations_voxelwise = var_innovations_voxelwise,
    na_mask = na_mask
  ))
}

.apply_ar_filter_to_matrix_cols <- function(M, ar_parameters_single_row, ar_order) {
  if (!is.numeric(M) || !is.matrix(M)) stop(".apply_ar_filter_to_matrix_cols: M must be a numeric matrix.")
  if (ar_order < 1) return(M)

  if (is.null(ar_parameters_single_row) || length(ar_parameters_single_row) != ar_order ||
      all(is.na(ar_parameters_single_row)) || all(ar_parameters_single_row == 0)){
      return(M)
  }

  n_timepoints_m <- nrow(M)
  if (n_timepoints_m <= ar_order) {
    warning(sprintf("Matrix to be filtered has %d rows, AR order is %d. Returning NA matrix.", n_timepoints_m, ar_order))
    return(matrix(NA_real_, nrow = n_timepoints_m, ncol = ncol(M)))
  }

  apply_ar_filter_matrix_cpp(M, ar_parameters_single_row)
}

.apply_ar_filter_voxelwise <- function(Y_data_matrix, voxel_ar_coeffs_matrix, ar_order) {
  if (!is.numeric(Y_data_matrix) || !is.matrix(Y_data_matrix)) stop(".apply_ar_filter_voxelwise: Y_data_matrix must be a numeric matrix.")
  if (ar_order < 1) return(Y_data_matrix)

  n_timepoints_y <- nrow(Y_data_matrix)
  if (n_timepoints_y <= ar_order) {
      warning(sprintf("Y_data has %d rows, AR order is %d. Returning NA matrix.", n_timepoints_y, ar_order))
      return(matrix(NA_real_, nrow = n_timepoints_y, ncol = ncol(Y_data_matrix)))
  }

  apply_ar_filter_voxelwise_cpp(Y_data_matrix, voxel_ar_coeffs_matrix)
}


#' Fit AR coefficients for a single voxel
#'
#' Helper function used by `ndx_ar2_whitening` to estimate AR coefficients and
#' innovation variance for one voxel.
#'
#' @param voxel_residuals Numeric vector of residuals for one voxel.
#' @param order Integer AR model order.
#' @param weights_vec Optional numeric vector of weights for weighted fitting.
#' @return Numeric vector containing the estimated AR coefficients followed by
#'   the innovation variance.
#' @keywords internal
.fit_ar_single_voxel <- function(voxel_residuals, order = 2L, weights_vec = NULL) {
  ar_fit <- NULL
  phi <- rep(0, order)
  var_pred <- NA_real_

  ar_fit_res <- tryCatch({
    if (stats::var(voxel_residuals, na.rm = TRUE) > .Machine$double.eps^0.5) {
      if (is.null(weights_vec)) {
        stats::ar.yw(voxel_residuals, aic = FALSE, order.max = order)
      } else {
        embed_mat <- stats::embed(voxel_residuals, order + 1)
        y_ar <- embed_mat[, 1]
        X_ar <- embed_mat[, -1, drop = FALSE]
        w_ar <- weights_vec[(order + 1):length(voxel_residuals)]
        if (sum(w_ar, na.rm = TRUE) > 0) {
          fit <- tryCatch(stats::lm.wfit(x = X_ar, y = y_ar, w = w_ar),
                          error = function(e) NULL)
          if (!is.null(fit) && length(fit$coefficients) == order) {
            list(ar = as.numeric(fit$coefficients),
                 var.pred = sum((fit$residuals^2) * w_ar) / sum(w_ar))
          } else {
            NULL
          }
        } else {
          NULL
        }
      }
    } else {
      NULL
    }
  }, error = function(e) {
    NULL
  })
  if (!is.null(ar_fit_res)) ar_fit <- ar_fit_res

  if (!is.null(ar_fit) && length(ar_fit$ar) == order) {
    phi_candidate <- ar_fit$ar
    roots <- tryCatch(polyroot(c(1, -phi_candidate)), error = function(e) NULL)
    if (!is.null(roots) && !any(is.na(roots)) && all(Mod(roots) > 1.00001)) {
      phi <- phi_candidate
      phi[abs(phi) < 0.01] <- 0
      var_pred <- ar_fit$var.pred
    }
  }

  c(phi, var_pred)
}

