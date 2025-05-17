# ND-X Utility Functions
# This file will contain helper functions for:
# - Loading fMRI data (NIfTI)
# - Basic preprocessing (demeaning, detrending if not part of Pass 0)
# - Other miscellaneous utilities

#' @import Rcpp neuroim2 oro.nifti psd rsvd
#' @importFrom stats convolve lm.fit
#' @importFrom oro.nifti origin reorient slice
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
  return(R2)
}
