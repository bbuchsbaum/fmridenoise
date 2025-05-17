#' Solve Ridge Regression for fMRI Data
#'
#' Solves the ridge regression problem for whitened fMRI data and a whitened design matrix.
#' For Sprint 1, this implements a basic isotropic ridge regression.
#'
#' @param Y_whitened A numeric matrix (timepoints x voxels) of whitened fMRI data.
#'   The first `order` rows (from AR pre-whitening) might be NA and should be handled (e.g. by removing them
#'   along with corresponding rows in X_whitened before solving).
#' @param X_whitened A numeric matrix (timepoints x regressors) of the whitened design matrix.
#'   Corresponding rows to NAs in Y_whitened should also be handled/removed.
#' @param lambda_ridge A single numeric value for the ridge penalty (lambda). For Sprint 1,
#'   this is a fixed value. GCV selection can be added later.
#' @param na_mask An optional logical vector indicating rows to remove from Y_whitened and X_whitened
#'   due to NA values (e.g., from AR filter initialization). If NULL (default), the function will
#'   attempt to identify and remove rows with any NAs in Y_whitened.
#'
#' @return A matrix (regressors x voxels) of estimated beta coefficients (betas_whitened).
#'
#' @details
#' The ridge regression solution is `beta = (X'X + lambda * I)^(-1) X'Y`.
#' This function assumes `Y_whitened` and `X_whitened` have had initial NA rows
#' (e.g., from AR pre-whitening) removed or properly handled if `na_mask` is provided.
#' If `na_mask` is not provided, any row in `Y_whitened` containing at least one NA will lead to that row (and the
#' corresponding row in `X_whitened`) being removed before model fitting.
#'
#' @export
ndx_solve_ridge <- function(Y_whitened, X_whitened, lambda_ridge, na_mask = NULL) {

  # Input validation
  if (!is.matrix(Y_whitened) || !is.numeric(Y_whitened)) {
    stop("Y_whitened must be a numeric matrix.")
  }
  if (!is.matrix(X_whitened) || !is.numeric(X_whitened)) {
    stop("X_whitened must be a numeric matrix.")
  }
  if (nrow(Y_whitened) != nrow(X_whitened)) {
    stop("Y_whitened and X_whitened must have the same number of rows (timepoints).")
  }
  if (!is.numeric(lambda_ridge) || length(lambda_ridge) != 1 || lambda_ridge < 0) {
    stop("lambda_ridge must be a single non-negative numeric value.")
  }

  n_timepoints_orig <- nrow(Y_whitened)
  n_voxels <- ncol(Y_whitened)
  n_regressors <- ncol(X_whitened)

  if (n_regressors == 0) {
    warning("X_whitened has 0 regressors. Returning empty beta matrix.")
    return(matrix(NA_real_, nrow = 0, ncol = n_voxels))
  }
  if (n_voxels == 0) {
    warning("Y_whitened has 0 voxels. Returning empty beta matrix.")
    return(matrix(NA_real_, nrow = n_regressors, ncol = 0))
  }

  # Handle NA rows from AR pre-whitening
  if (!is.null(na_mask)) {
    if (!is.logical(na_mask) || length(na_mask) != n_timepoints_orig) {
      stop("Provided na_mask must be a logical vector of length equal to nrows of Y_whitened.")
    }
    Y_clean <- Y_whitened[!na_mask, , drop = FALSE]
    X_clean <- X_whitened[!na_mask, , drop = FALSE]
  } else {
    # Default: remove rows with any NAs in Y_whitened (more robust if na_mask not passed)
    complete_y_rows <- stats::complete.cases(Y_whitened)
    if (!all(complete_y_rows)) {
        message(sprintf("Removing %d rows with NAs in Y_whitened prior to ridge regression.", sum(!complete_y_rows)))
    }
    Y_clean <- Y_whitened[complete_y_rows, , drop = FALSE]
    X_clean <- X_whitened[complete_y_rows, , drop = FALSE]
  }
  
  n_timepoints_clean <- nrow(X_clean)

  if (n_timepoints_clean == 0) {
    warning("No timepoints remaining after NA removal. Returning NA betas.")
    return(matrix(NA_real_, nrow = n_regressors, ncol = n_voxels))
  }
  if (n_timepoints_clean < n_regressors) {
    warning(sprintf("Number of timepoints after NA removal (%d) is less than number of regressors (%d). Ridge solution might be unstable or fail. Returning NA betas.", 
                    n_timepoints_clean, n_regressors))
    return(matrix(NA_real_, nrow = n_regressors, ncol = n_voxels))
  }

  # Ridge regression calculation: beta = (X'X + lambda * I)^-1 X'Y
  # More stable computation: solve((X'X + lambda*I) beta = X'Y)
  XtX <- crossprod(X_clean) # t(X_clean) %*% X_clean
  # Identity matrix of size n_regressors x n_regressors
  I <- diag(n_regressors) 
  
  # Ensure I has the same dimensions as XtX if n_regressors is 1
  if (n_regressors == 1 && !is.matrix(I)) {
    I <- matrix(I, 1, 1)
  }
  
  # LHS = (X'X + lambda * I)
  lhs <- XtX + lambda_ridge * I
  # RHS = X'Y
  XtY <- crossprod(X_clean, Y_clean) # t(X_clean) %*% Y_clean

  betas_whitened <- tryCatch({
    solve(lhs, XtY)
  }, error = function(e) {
    warning(paste("Solving ridge regression failed:", e$message, "Returning NA betas."))
    matrix(NA_real_, nrow = n_regressors, ncol = n_voxels)
  })
  
  return(betas_whitened)
}

#' Extract Task-Related Beta Coefficients
#'
#' Extracts beta coefficients corresponding to specified task regressors from the full
#' set of beta coefficients. For Sprint 1, this function primarily handles extraction.
#' The "unwhitening" aspect is a placeholder for future refinement.
#'
#' @param betas_whitened A matrix (total_regressors x voxels) of whitened beta coefficients,
#'   typically the output of `ndx_solve_ridge`.
#' @param X_whitened_colnames A character vector of column names for the whitened design matrix (`X_whitened`)
#'   that was used to estimate `betas_whitened`. This is used to identify task regressors by name.
#' @param task_regressor_names A character vector containing the names of the task regressors
#'   to be extracted.
#' @param ar_coeffs_global Optional. The global AR coefficients (vector of length `order`) that were used
#'   to whiten the design matrix `X`. If provided, a conceptual note about unwhitening can be considered.
#'   (Currently not used for actual unwhitening in Sprint 1).
#'
#' @return A matrix (num_task_regressors x voxels) of the extracted task-related beta coefficients.
#'   Row names will correspond to `task_regressor_names`.
#'   If no task regressors are found or an error occurs, returns NULL or an empty matrix with a warning.
#'
#' @details
#' Unwhitening beta coefficients is complex when different whitening strategies are applied
#' to Y and X (e.g., voxel-wise AR for Y, global AR for X). For Sprint 1, this function
#' focuses on the reliable extraction of beta coefficients based on their names.
#' A full unwhitening procedure that transforms betas back to the scale of an unwhitened model
#' (Y ~ X) may require further methodological development and is deferred.
#'
#' @export
ndx_extract_task_betas <- function(betas_whitened, X_whitened_colnames, task_regressor_names, ar_coeffs_global = NULL) {

  # Input validation
  if (!is.matrix(betas_whitened) || !is.numeric(betas_whitened)) {
    stop("betas_whitened must be a numeric matrix.")
  }
  if (nrow(betas_whitened) != length(X_whitened_colnames)) {
    stop("Number of rows in betas_whitened must match the length of X_whitened_colnames.")
  }
  if (!is.character(X_whitened_colnames)) {
    stop("X_whitened_colnames must be a character vector.")
  }
  if (!is.character(task_regressor_names) || length(task_regressor_names) == 0) {
    stop("task_regressor_names must be a non-empty character vector.")
  }

  n_total_regressors <- nrow(betas_whitened)
  n_voxels <- ncol(betas_whitened)

  if (n_total_regressors == 0) {
    warning("Input betas_whitened has 0 regressors. Returning NULL.")
    return(NULL)
  }
  
  # Match task_regressor_names to X_whitened_colnames
  # Ensure rownames of betas_whitened match X_whitened_colnames if they exist, otherwise use indices
  if (!is.null(rownames(betas_whitened)) && !identical(rownames(betas_whitened), X_whitened_colnames)) {
      # This case should ideally not happen if betas_whitened comes from a solve with X_whitened
      warning("Rownames of betas_whitened do not match X_whitened_colnames. Using X_whitened_colnames for indexing.")
  }
  
  # Find indices of task regressors
  task_indices <- match(task_regressor_names, X_whitened_colnames)
  
  if (any(is.na(task_indices))) {
    warning(sprintf("The following task regressors were not found in X_whitened_colnames: %s", 
                    paste(task_regressor_names[is.na(task_indices)], collapse=", ")))
    task_indices <- task_indices[!is.na(task_indices)] # Keep only found indices
  }

  if (length(task_indices) == 0) {
    warning("No specified task regressors were found in X_whitened_colnames. Returning NULL.")
    return(NULL)
  }
  
  extracted_betas <- betas_whitened[task_indices, , drop = FALSE]
  
  # Set rownames of extracted betas to the found task regressor names
  found_task_names <- X_whitened_colnames[task_indices]
  rownames(extracted_betas) <- found_task_names
  
  # Placeholder for unwhitening logic (Sprint 2+)
  if (!is.null(ar_coeffs_global)) {
    # message("Note: Unwhitening of betas is not implemented in Sprint 1. Betas are on the whitened scale.")
    # Conceptual: If X* = WX and Y* = WY (assume same W for simplicity for a moment)
    # Beta* from Y* = X*Beta* + e*
    # If W is invertible, Y = W^-1 Y*, X = W^-1 X*
    # How Beta* relates to Beta from Y = X Beta + e is complicated if W differs or is voxel-specific for Y.
  }
  
  return(extracted_betas)
} 