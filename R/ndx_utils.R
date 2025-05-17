# ND-X Utility Functions
# This file will contain helper functions for:
# - Loading fMRI data (NIfTI)
# - Basic preprocessing (demeaning, detrending if not part of Pass 0)
# - Other miscellaneous utilities

#' @import Rcpp neuroim2 oro.nifti psd rsvd
#' @importFrom stats convolve lm.fit
#' @importFrom oro.nifti origin reorient slice
#' @import matrixStats
NULL

# Placeholder for actual utility functions

#' Calculate OLS residuals using lm.fit
#'
#' @param Y Dependent variable matrix (timepoints x variables/voxels).
#' @param X Design matrix (timepoints x regressors).
#' @return Matrix of residuals (timepoints x variables/voxels).
#' @importFrom stats lm.fit
#' @keywords internal
#' @export
calculate_residuals_ols <- function(Y, X) {
  if (nrow(Y) != nrow(X)) {
    stop("Number of rows in Y and X must match for OLS.")
  }
  
  # If X has no columns (e.g. model with only intercept removed, or empty model)
  # then residuals are Y itself (or Y - mean(Y) if an intercept was implicitly fit and removed).
  # For an empty X, lm.fit might error or behave unexpectedly. Residuals are simply Y.
  if (ncol(X) == 0) {
    warning("Design matrix X has zero columns. Returning Y as residuals.")
    return(Y)
  }
  
  # Check for rank deficiency, lm.fit handles it by using a pseudo-inverse essentially
  qr_X <- qr(X)
  if (qr_X$rank < ncol(X)) {
    warning(paste("Design matrix is rank deficient. Rank =", qr_X$rank, "Columns =", ncol(X),
                  "OLS estimates will be non-unique for some regressors but residuals should be valid."))
  }

  # stats::lm.fit is efficient for multiple response variables (columns in Y)
  fit <- stats::lm.fit(X, Y)
  res <- fit$residuals
  # Remove any row or column names to ensure pure numeric matrix
  dimnames(res) <- NULL
  return(res)
}

#' Calculate Voxel-wise R-squared
#'
#' @param Y_observed Matrix of observed data (timepoints x voxels).
#' @param Y_residuals Matrix of residuals from a model (timepoints x voxels).
#' @return A numeric vector of R-squared values, one for each voxel.
#'   Names are stripped from the returned vector.
#' @details Valid observations for each voxel are selected using
#'   `!is.na(Y_observed[, j])`. Total Sum of Squares (TSS) and Residual Sum of
#'   Squares (RSS) are computed from these valid values only. If the calculated
#'   TSS is zero the returned R\eqn{^2} for that voxel is set to zero.
#' @importFrom matrixStats colVars
#' @keywords internal
#' @export
calculate_R2_voxelwise <- function(Y_observed, Y_residuals) {
  if (nrow(Y_observed) != nrow(Y_residuals) || ncol(Y_observed) != ncol(Y_residuals)) {
    stop("Dimensions of Y_observed and Y_residuals must match.")
  }

  n_vox <- ncol(Y_observed)
  R2 <- numeric(n_vox)
  for (j in seq_len(n_vox)) {
    mask <- !is.na(Y_observed[, j])
    if (sum(mask) == 0) {
      R2[j] <- 0
      next
    }

    y_obs <- Y_observed[mask, j]
    y_res <- Y_residuals[mask, j]
    TSS <- sum((y_obs - mean(y_obs))^2)
    RSS <- sum(y_res^2, na.rm = TRUE)

    if (abs(TSS) < 1e-9) {
      R2[j] <- 0
    } else {
      r2_val <- 1 - (RSS / TSS)
      if (!is.finite(r2_val) || r2_val < 0) r2_val <- 0
      R2[j] <- r2_val
    }
  }

  R2 <- unname(R2)
  R2[R2 < 0] <- 0 
  return(R2)
}

#' Calculate Denoising Efficacy Score (DES)
#'
#' DES = 1 - (Variance of current residuals / Variance of baseline residuals)
#'
#' @param current_residuals_unwhitened Numeric matrix or vector of residuals from the current denoising model.
#' @param VAR_BASELINE_FOR_DES Numeric scalar, the variance of residuals from a baseline model 
#'   (e.g., task-only GLM, or residuals from Pass 0 of ndx_initial_glm).
#' @return Numeric scalar, the Denoising Efficacy Score. Returns NA if inputs are invalid.
#' @importFrom stats var
#' @export
calculate_DES <- function(current_residuals_unwhitened, VAR_BASELINE_FOR_DES) {
  if (is.null(current_residuals_unwhitened) || is.null(VAR_BASELINE_FOR_DES)) {
    warning("calculate_DES: Inputs cannot be NULL.")
    return(NA_real_)
  }
  if (!is.numeric(current_residuals_unwhitened) || !is.numeric(VAR_BASELINE_FOR_DES)) {
    warning("calculate_DES: Inputs must be numeric.")
    return(NA_real_)
  }
  if (length(VAR_BASELINE_FOR_DES) != 1 || !is.finite(VAR_BASELINE_FOR_DES)){
    warning("calculate_DES: VAR_BASELINE_FOR_DES must be a single finite number.")
    return(NA_real_)
  }
  if (VAR_BASELINE_FOR_DES <= 1e-9) { # Avoid division by zero or near-zero
    warning("calculate_DES: VAR_BASELINE_FOR_DES is too close to zero. DES calculation may be unstable or meaningless.")
    # Depending on convention, could return NA, 0, or 1 if current_residuals also have ~0 variance.
    # If baseline variance is 0, any non-zero current residual variance would yield DES < 0 (formally -Inf).
    # If both are 0, DES is undefined (0/0) or could be 1 (perfect explanation by baseline).
    # For now, return NA as it's likely an issue with baseline variance calculation.
    return(NA_real_)
  }

  var_current_residuals <- stats::var(as.vector(current_residuals_unwhitened), na.rm = TRUE)
  if (!is.finite(var_current_residuals)){
      warning("calculate_DES: Variance of current_residuals_unwhitened is not finite.")
      return(NA_real_)
  }

  des <- 1 - (var_current_residuals / VAR_BASELINE_FOR_DES)
  return(des)
}

#' Orthogonalize a Target Matrix Against a Basis Matrix
#'
#' This function orthogonalizes each column of a `target_matrix` with respect to
#' all columns in a `basis_matrix`. The orthogonalization is performed by
#' successively regressing each target column against the `basis_matrix` and
#' taking the residuals. This is a form of Gram-Schmidt orthogonalization.
#'
#' @param target_matrix A numeric matrix whose columns will be orthogonalized.
#' @param basis_matrix A numeric matrix representing the basis against which
#'   the columns of `target_matrix` will be orthogonalized. Must have the same
#'   number of rows as `target_matrix`.
#' @param tol Numeric, tolerance for checking if a basis vector is zero or if a
#'  target vector becomes zero after orthogonalization. Columns in `basis_matrix`
#'  with sum of squares less than `tol` will be excluded. Target columns that become
#'  zero (sum of squares < `tol`) after orthogonalization will be returned as zeros.
#'  Default is `1e-9`.
#'
#' @return A matrix with the same dimensions as `target_matrix`, where each
#'   column is the orthogonalized version of the corresponding column in the
#'   input `target_matrix`. If `basis_matrix` is NULL, has zero columns, or all its
#'   columns have near-zero variance, the original `target_matrix` is returned.
#'   If `target_matrix` is NULL or has zero columns, it is returned as is.
#'
#' @examples
#' \dontrun{
#' basis <- matrix(rnorm(10*2), 10, 2)
#' basis[,2] <- basis[,1] * 2 + rnorm(10, 0, 0.1) # Make somewhat collinear
#' target <- matrix(rnorm(10*3), 10, 3)
#' target[,1] <- basis[,1] * 0.5 + rnorm(10, 0.5)
#' target[,2] <- basis[,2] * 0.8 + rnorm(10, 0.2)
#'
#' orthogonalized_target <- ndx_orthogonalize_matrix_against_basis(target, basis)
#'
#' # Check orthogonality (should be close to zero)
#' # t(orthogonalized_target) %*% basis
#'
#' # Example with a zero column in basis
#' basis_with_zero <- cbind(basis, rep(0, 10))
#' ortho_target_zero_basis <- ndx_orthogonalize_matrix_against_basis(target, basis_with_zero)
#'
#' # Example with basis_matrix being NULL
#' ortho_target_null_basis <- ndx_orthogonalize_matrix_against_basis(target, NULL)
#' print(all.equal(target, ortho_target_null_basis))
#' }
#' @export
ndx_orthogonalize_matrix_against_basis <- function(target_matrix, basis_matrix, tol = 1e-9) {

  if (is.null(target_matrix) || ncol(target_matrix) == 0) {
    return(target_matrix)
  }
  if (!is.matrix(target_matrix) || !is.numeric(target_matrix)){
    stop("`target_matrix` must be a numeric matrix.")
  }

  if (is.null(basis_matrix) || ncol(basis_matrix) == 0) {
    return(target_matrix)
  }
  if (!is.matrix(basis_matrix) || !is.numeric(basis_matrix)){
    stop("`basis_matrix` must be a numeric matrix.")
  }

  if (nrow(target_matrix) != nrow(basis_matrix)) {
    stop("`target_matrix` and `basis_matrix` must have the same number of rows.")
  }

  # Filter out near-zero variance columns from basis_matrix
  valid_basis_cols_ss <- colSums(basis_matrix^2, na.rm = TRUE)
  valid_basis_indices <- which(valid_basis_cols_ss >= tol)

  if (length(valid_basis_indices) == 0) {
    warning("All columns in `basis_matrix` have near-zero variance. Returning original `target_matrix`.")
    return(target_matrix)
  }
  effective_basis_matrix <- basis_matrix[, valid_basis_indices, drop = FALSE]

  orthogonalized_matrix <- matrix(0.0, nrow = nrow(target_matrix), ncol = ncol(target_matrix))
  colnames(orthogonalized_matrix) <- colnames(target_matrix)
  rownames(orthogonalized_matrix) <- rownames(target_matrix)

  for (i in 1:ncol(target_matrix)) {
    y_target_col <- target_matrix[, i, drop = FALSE]
    
    # Check if target column itself is near zero
    if (sum(y_target_col^2, na.rm = TRUE) < tol) {
      orthogonalized_matrix[, i] <- y_target_col # Already zero (or near zero)
      next
    }

    # Using lm.fit for efficiency to get residuals
    # y_target_col ~ effective_basis_matrix
    # Residuals are y_target_col - X_basis * ( (t(X_basis) %*% X_basis)^-1 %*% t(X_basis) %*% y_target_col )
    # No intercept is added by lm.fit by default, which is what we want here.
    fit <- tryCatch({
      stats::lm.fit(x = effective_basis_matrix, y = y_target_col)
    }, error = function(e) {
      warning(sprintf("lm.fit failed for target column %d during orthogonalization: %s. Column set to original.", i, e$message))
      return(NULL) # Indicate failure
    })

    if (!is.null(fit) && !is.null(fit$residuals)) {
      residuals_ortho <- fit$residuals
      if (sum(residuals_ortho^2, na.rm = TRUE) < tol) {
        # If residual is near zero, it means the target column was fully explained by the basis
        # We keep it as zero vector (already initialized in orthogonalized_matrix)
        # orthogonalized_matrix[, i] <- rep(0.0, nrow(target_matrix))
      } else {
        orthogonalized_matrix[, i] <- residuals_ortho
      }
    } else {
      # If lm.fit failed, or did not produce residuals, keep original target column
      # A warning would have been issued by the tryCatch error handler
      orthogonalized_matrix[, i] <- y_target_col
    }
  }

  return(orthogonalized_matrix)
}

#' Merge Two Lists with Defaults
#'
#' Wrapper around `utils::modifyList` that also handles `NULL` user lists.
#'
#' @param defaults A list of default values.
#' @param user User-supplied options.
#' @return Combined list.
#' @keywords internal
#' @export
merge_lists <- function(defaults, user) {
  utils::modifyList(defaults, if (is.null(user)) list() else user)
}

#' Null-Coalescing Operator
#'
#' Return `b` when `a` is `NULL`, otherwise `a`.
#'
#' @param a First value.
#' @param b Fallback value.
#' @keywords internal
#' @export
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}

