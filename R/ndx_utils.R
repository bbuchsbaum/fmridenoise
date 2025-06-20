# ND-X Utility Functions
# This file will contain helper functions for:
# - Loading fMRI data (NIfTI)
# - Basic preprocessing (demeaning, detrending if not part of Pass 0)
# - Other miscellaneous utilities

#' @useDynLib ndx, .registration = TRUE
#' @import Rcpp neuroim2 oro.nifti psd rsvd
#' @importFrom stats convolve lm.fit
#' @importFrom oro.nifti origin reorient slice
#' @import matrixStats
NULL


# Placeholder for actual utility functions

#' Return `NA` if Value is `NULL`
#'
#' Simple helper that returns `NA` when the input is `NULL`, otherwise the
#' input value unchanged.
#'
#' @param x An object that might be `NULL`.
#' @return `NA` if `x` is `NULL`, otherwise `x`.
#' @export
ndx_val_or_na <- function(x) {
  if (is.null(x)) NA else x
}

#' Calculate OLS residuals using lm.fit
#'
#' @param Y Dependent variable matrix (timepoints x variables/voxels).
#' @param X Design matrix (timepoints x regressors).
#' @return Matrix of residuals (timepoints x variables/voxels).
#' @importFrom stats lm.fit
#' @export
calculate_residuals_ols <- function(Y, X) {
  # Convert inputs to numeric matrices if possible
  if (!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if (!is.numeric(Y) || !is.numeric(X)) {
    stop("`Y` and `X` must be numeric matrices.")
  }

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
#' @export
calculate_R2_voxelwise <- function(Y_observed, Y_residuals) {
  if (nrow(Y_observed) != nrow(Y_residuals) || ncol(Y_observed) != ncol(Y_residuals)) {
    stop("Dimensions of Y_observed and Y_residuals must match.")
  }

  mask <- !is.na(Y_observed)

  valid_counts <- matrixStats::colSums2(mask)

  obs_clean <- Y_observed
  obs_clean[!mask] <- 0

  means <- matrixStats::colSums2(obs_clean) / ifelse(valid_counts > 0, valid_counts, 1)

  centered <- (obs_clean - matrix(means, nrow = nrow(Y_observed), ncol = ncol(Y_observed), byrow = TRUE)) * mask
  TSS <- matrixStats::colSums2(centered^2)

  res_clean <- Y_residuals
  res_clean[is.na(res_clean)] <- 0
  RSS <- matrixStats::colSums2((res_clean^2) * mask)

  R2 <- ifelse(valid_counts == 0 | TSS < 1e-9, 0, 1 - RSS / TSS)
  R2[!is.finite(R2) | R2 < 0] <- 0

  unname(R2)
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
  if (!is.numeric(current_residuals_unwhitened)) {
    current_residuals_unwhitened <- as.numeric(current_residuals_unwhitened)
    if (anyNA(current_residuals_unwhitened)) {
      warning("calculate_DES: `current_residuals_unwhitened` contains non-numeric values.")
      return(NA_real_)
    }
  }
  if (!is.numeric(VAR_BASELINE_FOR_DES)) {
    warning("calculate_DES: `VAR_BASELINE_FOR_DES` must be numeric.")
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
#' all columns in a `basis_matrix`. The orthogonalization is performed using a
#' projection onto the orthonormal basis of `basis_matrix` computed via QR
#' decomposition. Formally the result is
#' `target_matrix - Q %*% crossprod(Q, target_matrix)` where `Q` contains the
#' orthonormal columns spanning `basis_matrix`.
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
  if (!is.matrix(target_matrix)) {
    target_matrix <- as.matrix(target_matrix)
  }
  if (!is.numeric(target_matrix)) {
    stop("`target_matrix` must be a numeric matrix.")
  }

  if (is.null(basis_matrix) || ncol(basis_matrix) == 0) {
    return(target_matrix)
  }
  if (!is.matrix(basis_matrix)) {
    basis_matrix <- as.matrix(basis_matrix)
  }
  if (!is.numeric(basis_matrix)) {
    stop("`basis_matrix` must be a numeric matrix.")
  }

  if (nrow(target_matrix) != nrow(basis_matrix)) {
    stop("`target_matrix` and `basis_matrix` must have the same number of rows.")
  }

  # Filter out near-zero variance columns from basis_matrix
  valid_basis_cols_ss <- colSums(basis_matrix^2, na.rm = TRUE)
  valid_basis_indices <- which(valid_basis_cols_ss >= tol)

  if (length(valid_basis_indices) == 0) {
    return(target_matrix)
  }

  B <- basis_matrix[, valid_basis_indices, drop = FALSE]
  qr_basis <- qr(B)
  Q <- qr.Q(qr_basis)

  result <- target_matrix - Q %*% crossprod(Q, target_matrix)
  colnames(result) <- colnames(target_matrix)
  rownames(result) <- rownames(target_matrix)
  result
}

#' Merge Two Lists with Defaults
#'
#' Wrapper around `utils::modifyList` that also handles `NULL` user lists.
#'
#' @param defaults A list of default values.
#' @param user User-supplied options.
#' @return Combined list.
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
#' @return `a` if not `NULL`, otherwise `b`.
#' @keywords internal
#' @export
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}

#' Calculate Beta Stability Across Runs
#'
#' Computes average pairwise correlations of task beta coefficients across runs
#' for each pass of the ND-X workflow.
#'
#' @param betas_per_pass A list where each element corresponds to a pass and
#'   contains a list of beta matrices (one per run).
#' @return Numeric vector of average cross-run correlations for each pass. If a
#'   pass contains fewer than two runs, the returned value is `NA`.
#' @export
calculate_beta_stability <- function(betas_per_pass) {
  if (!is.list(betas_per_pass) || length(betas_per_pass) == 0) {
    stop("betas_per_pass must be a non-empty list")
  }
  stability <- numeric(length(betas_per_pass))
  for (i in seq_along(betas_per_pass)) {
    run_betas <- betas_per_pass[[i]]
    if (!is.list(run_betas) || length(run_betas) < 2) {
      stability[i] <- NA_real_
      next
    }
    vecs <- lapply(run_betas, function(b) as.vector(b))
    cor_vals <- utils::combn(seq_along(vecs), 2, function(idx) {
      stats::cor(vecs[[idx[1]]], vecs[[idx[2]]], use = "pairwise.complete.obs")
    })
    stability[i] <- mean(cor_vals, na.rm = TRUE)
  }
  stability
}

#' Compute Ljung-Box p-values for Residuals
#'
#' Applies the Ljung-Box test to each residual time series and returns the
#' p-values.
#'
#' @param residuals Matrix of residuals (timepoints x voxels).
#' @param lag Number of lags to use in the Ljung-Box test. Default 5.
#' @return Numeric vector of p-values, one per column of `residuals`.
#' @export
compute_ljung_box_pvalues <- function(residuals, lag = 5L) {
  if (!is.matrix(residuals) || !is.numeric(residuals)) {
    stop("residuals must be a numeric matrix")
  }
  apply(residuals, 2, function(ts) {
    if (all(is.na(ts))) return(NA_real_)
    # Ensure there are enough non-NA observations for the test
    if (sum(!is.na(ts)) <= lag) return(NA_real_)
    stats::Box.test(ts, lag = lag, type = "Ljung-Box")$p.value
  })
}

#' Compute Ljung-Box p-value for Residual Whiteness
#'
#' Given a matrix or vector of residuals, this helper collapses the data across
#' voxels (columns) by averaging and then applies `stats::Box.test` with the
#' Ljung-Box option to quantify remaining autocorrelation.
#'
#' @param residuals Numeric matrix (timepoints x voxels) or vector of residuals.
#' @param lag Integer lag parameter passed to `stats::Box.test`. Default 5.
#' @return Numeric p-value from the Ljung-Box test or `NA` if it cannot be
#'   computed.
#' @export
ndx_ljung_box_pval <- function(residuals, lag = 5L) {
  if (is.null(residuals)) return(NA_real_)
  if (is.vector(residuals)) {
    vec <- as.numeric(residuals)
  } else if (is.matrix(residuals)) {
    vec <- rowMeans(residuals, na.rm = TRUE)
  } else {
    # stop("residuals must be a numeric vector or matrix") # Or handle differently
    return(NA_real_)
  }

  vec <- vec[is.finite(vec)] # Remove NAs and Infs
  if (length(vec) <= lag) return(NA_real_) # Check after NA removal

  # Check for zero variance, which can cause Box.test to fail
  if (stats::var(vec, na.rm = TRUE) < .Machine$double.eps * 100) return(NA_real_)

  tryCatch({
    stats::Box.test(vec, lag = lag, type = "Ljung-Box")$p.value
  }, error = function(e) {
    warning(paste("Ljung-Box test failed:", e$message))
    NA_real_
  })
}

#' Calculate Annihilation Mode Verdict Statistics
#'
#' This utility computes the variance ratio between ND-X unique components and
#' GLMdenoise components and returns the categorical verdict used throughout the
#' ND-X diagnostics.
#'
#' @param workflow_output List produced by `NDX_Process_Subject` containing the
#'   matrices `gdlite_pcs`, `rpca_orthogonalized`, and
#'   `spectral_orthogonalized`.
#' @return A list with elements `var_ratio` and `verdict`. If necessary inputs
#'   are missing or variance of GLMdenoise components is zero, returns NA values.
#' @export
ndx_annihilation_verdict_stats <- function(workflow_output) {
  # Ensure all required components are present
  required_gdlite <- !is.null(workflow_output$gdlite_pcs)
  required_ndx_unique <- !is.null(workflow_output$rpca_orthogonalized) ||
                         !is.null(workflow_output$spectral_orthogonalized)

  if (!required_gdlite || !required_ndx_unique) {
    return(list(var_ratio = NA_real_, verdict = NA_character_))
  }

  var_gd <- sum(workflow_output$gdlite_pcs^2, na.rm = TRUE)

  # If GLMdenoise components have zero variance, ratio is undefined or infinite.
  # Return NA as it's not a meaningful comparison.
  if (var_gd <= .Machine$double.eps) {
    return(list(var_ratio = NA_real_, verdict = NA_character_))
  }
  
  var_ndx <- 0
  if (!is.null(workflow_output$rpca_orthogonalized)) {
    var_ndx <- var_ndx + sum(workflow_output$rpca_orthogonalized^2, na.rm = TRUE)
  }
  if (!is.null(workflow_output$spectral_orthogonalized)) {
    var_ndx <- var_ndx + sum(workflow_output$spectral_orthogonalized^2, na.rm = TRUE)
  }

  ratio <- var_ndx / var_gd
  verdict <- if (ratio < 0.1) {
    "Tie"
  } else if (ratio < 0.5) {
    "Win"
  } else if (ratio < 1.0) {
    "Decisive Win"
  } else { # ratio >= 1.0
    "Annihilation"
  }

  list(var_ratio = ratio, verdict = verdict)
}

