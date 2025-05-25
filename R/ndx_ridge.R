#' Compute Projection Matrices for Anisotropic Ridge
#'
#' Given orthonormal basis matrices describing different nuisance subspaces,
#' this helper constructs the projection matrices used by anisotropic ridge
#' regression. Any NULL or empty basis results in a NULL projector.
#' This function is part of the modular anisotropic ridge API.
#'
#' @param U_GD Matrix of GLMdenoise PCs (n_regressors x k1). Optional.
#' @param U_Unique Matrix of nuisance components unique to ND-X (n_regressors x k2). Optional.
#' @param U_Noise Matrix of nuisance components when not using Annihilation Mode (n_regressors x k3). Optional.
#' @param n_regressors Integer number of regressors in the design.
#' @return A list with projection matrices `P_GD`, `P_Unique`, `P_Noise`, and `P_Signal` (the complement).
#' @export
ndx_compute_projection_matrices <- function(U_GD = NULL, U_Unique = NULL,
                                            U_Noise = NULL, n_regressors) {
  I_reg <- diag(1, n_regressors)
  make_proj <- function(U) {
    if (is.null(U) || ncol(U) == 0) return(NULL)
    Q <- qr.Q(qr(U))
    Q %*% t(Q)
  }
  P_GD <- make_proj(U_GD)
  P_Unique <- make_proj(U_Unique)
  P_Noise <- make_proj(U_Noise)

  basis_list <- list(U_GD, U_Unique, U_Noise)
  non_null <- Filter(function(x) !is.null(x) && ncol(x) > 0, basis_list)
  if (length(non_null) > 0) {
    Q_all <- qr.Q(qr(do.call(cbind, non_null)))
    P_combined <- Q_all %*% t(Q_all)
  } else {
    P_combined <- matrix(0, n_regressors, n_regressors)
  }
  P_Signal <- I_reg - P_combined

  list(P_GD = P_GD, P_Unique = P_Unique, P_Noise = P_Noise, P_Signal = P_Signal)
}

#' Tune Parallel Lambda via Generalized Cross Validation
#'
#' Performs a simple GCV search over a grid of lambda_parallel values. The
#' lambda for the orthogonal complement (`lambda_perp`) is set as
#' `lambda_ratio * lambda_parallel`.
#' This helper is part of the modular anisotropic ridge API.
#'
#' @param Y_whitened Numeric matrix of whitened responses.
#' @param X_whitened Numeric matrix of whitened design.
#' @param P_Noise Projection matrix onto the noise subspace.
#' @param lambda_grid Numeric vector of candidate lambda_parallel values.
#' @param lambda_ratio Numeric multiplier for lambda_perp.
#' @return Numeric scalar best lambda_parallel from the grid.
#' @export
ndx_gcv_tune_lambda_parallel <- function(Y_whitened, X_whitened, P_Noise,
                                         lambda_grid = 10^seq(-2, 2, length.out = 5),
                                         lambda_ratio = 0.05) {
  n <- nrow(X_whitened)
  I_reg <- diag(1, ncol(X_whitened))
  if (is.null(P_Noise)) {
    diag_noise <- rep(0, ncol(X_whitened))
    P_Signal <- I_reg
  } else {
    diag_noise <- diag(P_Noise)
    P_Signal <- I_reg - P_Noise
  }
  best_lambda <- lambda_grid[1]
  best_gcv <- Inf
  XtX <- crossprod(X_whitened)
  XtY <- crossprod(X_whitened, Y_whitened)
  base_diag <- diag(XtX)
  for (lam in lambda_grid) {
    K_diag <- lam * diag_noise + (lam * lambda_ratio) * diag(P_Signal)
    lhs <- XtX
    diag(lhs) <- base_diag + pmax(K_diag, .Machine$double.eps)
    chol_decomp <- tryCatch(chol(lhs), error = function(e) NULL)
    if (is.null(chol_decomp)) next
    beta_hat <- chol2inv(chol_decomp) %*% XtY
    residuals <- Y_whitened - X_whitened %*% beta_hat
    sse <- sum(residuals^2)
    S_mat <- X_whitened %*% chol2inv(chol_decomp) %*% t(X_whitened)
    edf <- sum(diag(S_mat))
    if ((n - edf) <= .Machine$double.eps) {
      gcv <- Inf
    } else {
      gcv <- sse / (n - edf)^2
    }
    if (is.finite(gcv) && gcv < best_gcv) {
      best_gcv <- gcv
      best_lambda <- lam
    }
  }
  best_lambda
}

#' Tune Ridge Penalty via Leave-One-Out Cross Validation
#'
#' Computes the optimal isotropic ridge penalty using the efficient
#' leave-one-out (LOO) formula. The cross-validation error is
#' calculated for each candidate value in `lambda_grid` and the
#' lambda yielding the smallest mean squared LOO residual is returned.
#'
#' @param Y Numeric matrix of responses (timepoints x variables).
#' @param X Numeric design matrix (timepoints x regressors).
#' @param lambda_grid Numeric vector of candidate ridge penalties.
#' @return The lambda value from `lambda_grid` with minimal LOO error.
#' @export
ndx_loocv_tune_lambda_ridge <- function(Y, X, lambda_grid = 10^seq(-2, 2, length.out = 5)) {
  if (!is.matrix(Y) || !is.numeric(Y)) stop("Y must be a numeric matrix.")
  if (!is.matrix(X) || !is.numeric(X)) stop("X must be a numeric matrix.")
  if (nrow(Y) != nrow(X)) stop("Y and X must have the same number of rows.")
  if (!is.numeric(lambda_grid) || length(lambda_grid) == 0) stop("lambda_grid must be numeric.")

  m <- ncol(X)
  XtX <- crossprod(X)
  XtY <- crossprod(X, Y)

  cv_errors <- numeric(length(lambda_grid))

  for (i in seq_along(lambda_grid)) {
    lam <- lambda_grid[i]
    cho <- chol(XtX + diag(lam, m))
    beta_all <- backsolve(cho, backsolve(cho, XtY, transpose = TRUE))
    inv_cho <- chol2inv(cho)
    hat_diag <- rowSums((X %*% inv_cho) * X)
    loo_res <- (Y - X %*% beta_all) / (1 - hat_diag)
    cv_errors[i] <- mean(loo_res^2)
  }

  lambda_grid[which.min(cv_errors)]
}

#' Update Lambda Aggressiveness Based on Rho
#'
#' Simple heuristic to adjust lambda_parallel depending on how much residual
#' variance projects onto the noise subspace.
#'
#' Part of the modular anisotropic ridge API.
#'
#' @param lambda_parallel Current lambda_parallel value.
#' @param rho_noise_projection Proportion of residual variance projected onto the noise subspace.
#' @param high_threshold Increase lambda if rho > this value. Default 0.10.
#' @param low_threshold Decrease lambda if rho < this value. Default 0.03.
#' @param increase_factor Multiplicative increase factor. Default 1.25.
#' @param decrease_factor Multiplicative decrease factor. Default 0.8.
#' @return Updated lambda_parallel.
#' @export
ndx_update_lambda_aggressiveness <- function(lambda_parallel, rho_noise_projection,
                                             high_threshold = 0.10, low_threshold = 0.03,
                                             increase_factor = 1.25, decrease_factor = 0.8) {
  if (is.na(rho_noise_projection)) return(lambda_parallel)
  if (rho_noise_projection > high_threshold) {
    lambda_parallel <- lambda_parallel * increase_factor
  } else if (rho_noise_projection < low_threshold) {
    lambda_parallel <- lambda_parallel * decrease_factor
  }
  lambda_parallel
}

#' Estimate Residual Variance after Whitening
#'
#' Computes a simple estimate of the residual variance from a whitened
#' design matrix and response matrix. This is used for scaling the lambda
#' grid during GCV search in anisotropic ridge regression.
#' Part of the modular anisotropic ridge API.
#'
#' @param Y_whitened Numeric matrix of whitened responses.
#' @param X_whitened Numeric matrix of whitened design.
#' @param na_mask Optional logical vector identifying rows to exclude
#'   (typically the initial rows lost during AR whitening).
#' @return Numeric scalar variance estimate. Returns `NA` if the estimate
#'   cannot be computed.
#' @export
ndx_estimate_res_var_whitened <- function(Y_whitened, X_whitened, na_mask = NULL) {
  if (!is.null(na_mask)) {
    Y_eff <- Y_whitened[!na_mask, , drop = FALSE]
    X_eff <- X_whitened[!na_mask, , drop = FALSE]
  } else {
    Y_eff <- Y_whitened
    X_eff <- X_whitened
  }

  if (nrow(X_eff) == 0 || ncol(X_eff) == 0) return(NA_real_)

  fit <- tryCatch(stats::lm.fit(x = X_eff, y = Y_eff), error = function(e) NULL)
  if (is.null(fit) || is.null(fit$residuals)) return(NA_real_)

  stats::var(as.vector(fit$residuals), na.rm = TRUE)
}

#' Validate Inputs for Ridge Solver
#'
#' Performs the set of sanity checks originally embedded in
#' `ndx_solve_anisotropic_ridge`.
#'
#' @keywords internal
validate_ridge_inputs <- function(Y_whitened, X_whitened, K_penalty_diag = NULL,
                                  na_mask = NULL, projection_mats = NULL,
                                  lambda_values = NULL, weights = NULL,
                                  K_penalty_mat = NULL,
                                  use_penalty_matrix = FALSE) {
  if (is.null(Y_whitened) || is.null(X_whitened)) {
    warning("Y_whitened and X_whitened must be provided.")
    return(FALSE)
  }
  if (!is.matrix(Y_whitened) || !is.numeric(Y_whitened)) {
    warning("Y_whitened must be a numeric matrix.")
    return(FALSE)
  }
  if (!is.matrix(X_whitened) || !is.numeric(X_whitened)) {
    warning("X_whitened must be a numeric matrix.")
    return(FALSE)
  }
  if (nrow(Y_whitened) != nrow(X_whitened)) {
    warning("Y_whitened and X_whitened must have the same number of rows (timepoints).")
    return(FALSE)
  }

  if (!is.null(weights)) {
    if (is.matrix(weights)) {
      if (!is.numeric(weights) || nrow(weights) != nrow(Y_whitened) ||
          ncol(weights) != ncol(Y_whitened)) {
        warning("weights matrix must match dimensions of Y_whitened")
        return(FALSE)
      }
      if (any(!is.finite(as.vector(weights))) || any(as.vector(weights) < 0)) {
        warning("weights must contain non-negative finite numbers")
        return(FALSE)
      }
    } else if (is.numeric(weights) && is.vector(weights)) {
      if (length(weights) != nrow(Y_whitened)) {
        warning("weights vector length must equal number of rows in Y_whitened")
        return(FALSE)
      }
      if (any(!is.finite(weights)) || any(weights < 0)) {
        warning("weights must contain non-negative finite numbers")
        return(FALSE)
      }
    } else {
      warning("weights must be a numeric vector or matrix")
      return(FALSE)
    }
  }

  if (use_penalty_matrix) {
    if (!is.null(K_penalty_mat)) {
      if (!is.matrix(K_penalty_mat) ||
          nrow(K_penalty_mat) != ncol(X_whitened) ||
          ncol(K_penalty_mat) != ncol(X_whitened)) {
        warning("K_penalty_mat must be a square matrix with dimension equal to ncol(X_whitened).")
        return(FALSE)
      }
    } else {
      if (is.null(projection_mats) || length(lambda_values) == 0) {
        warning("Either K_penalty_mat or projection_mats and lambda_values must be provided when use_penalty_matrix = TRUE.")
        return(FALSE)
      }
      if (!is.list(projection_mats) || is.null(projection_mats$P_Signal)) {
        warning("projection_mats must contain at least P_Signal.")
        return(FALSE)
      }
    }
  } else {
    if (!is.null(K_penalty_diag)) {
      if (!is.numeric(K_penalty_diag) || length(K_penalty_diag) != ncol(X_whitened)) {
        warning("K_penalty_diag must be a numeric vector with length equal to ncol(X_whitened).")
        return(FALSE)
      }
      if (any(K_penalty_diag < 0)) {
        warning("All elements of K_penalty_diag must be non-negative.")
        return(FALSE)
      }
    } else {
      if (is.null(projection_mats) || length(lambda_values) == 0) {
        warning("Either K_penalty_diag or projection_mats and lambda_values must be provided.")
        return(FALSE)
      }
      if (!is.list(projection_mats) || is.null(projection_mats$P_Signal)) {
        warning("projection_mats must contain at least P_Signal.")
        return(FALSE)
      }
    }
  }
  TRUE
}

#' Apply NA Mask and Weights
#'
#' Removes rows indicated by `na_mask` and subsets `weights` accordingly.
#'
#' @keywords internal
apply_na_mask_and_weights <- function(Y, X, na_mask = NULL, weights = NULL) {
  if (!is.null(na_mask)) {
    if (!is.logical(na_mask) || length(na_mask) != nrow(Y)) {
      warning("na_mask must be a logical vector with length equal to nrow(Y_whitened). Using all timepoints.")
      Y_eff <- Y
      X_eff <- X
      W_eff <- if (is.null(weights)) NULL else weights
    } else {
      if (all(na_mask)) {
        warning("All timepoints are masked by na_mask. Cannot solve ridge regression.")
        return(NULL)
      }
      Y_eff <- Y[!na_mask, , drop = FALSE]
      X_eff <- X[!na_mask, , drop = FALSE]
      if (!is.null(weights)) {
        if (is.matrix(weights)) {
          W_eff <- weights[!na_mask, , drop = FALSE]
        } else {
          W_eff <- weights[!na_mask]
        }
      } else {
        W_eff <- NULL
      }
    }
  } else {
    Y_eff <- Y
    X_eff <- X
    W_eff <- if (is.null(weights)) NULL else weights
  }

  if (nrow(X_eff) == 0) {
    warning("No timepoints remaining after NA removal (or initial input had 0 timepoints). Cannot solve.")
    return(NULL)
  }

  if (nrow(X_eff) < ncol(X_eff)) {
    warning(sprintf("Number of timepoints after NA removal (%d) is less than number of regressors (%d). Ridge solution might be unstable or not meaningful.",
                    nrow(X_eff), ncol(X_eff)))
  }

  list(Y = Y_eff, X = X_eff, W = W_eff)
}

#' Construct Ridge Penalty
#'
#' Builds the penalty matrix or diagonal from projection matrices and lambda values.
#'
#' @keywords internal
construct_penalty <- function(X_eff, Y_eff, K_penalty_diag = NULL,
                              projection_mats = NULL, lambda_values = NULL,
                              gcv_lambda = FALSE, res_var_scale = 1.0,
                              lambda_grid = 10^seq(-2, 2, length.out = 5),
                              K_penalty_mat = NULL, use_penalty_matrix = FALSE) {

  if (use_penalty_matrix) {
    if (is.null(K_penalty_mat)) {
      lambda_parallel <- lambda_values$lambda_parallel %||% 1.0
      lambda_signal   <- lambda_values$lambda_perp_signal %||% (0.05 * lambda_parallel)
      lambda_gd       <- lambda_values$lambda_gd %||% lambda_parallel
      lambda_unique   <- lambda_values$lambda_unique %||% lambda_parallel

      if (gcv_lambda) {
        if (!is.null(projection_mats$P_Noise)) {
          lambda_ratio <- lambda_signal / lambda_parallel
          lambda_parallel <- ndx_gcv_tune_lambda_parallel(Y_eff, X_eff,
                                                          projection_mats$P_Noise,
                                                          lambda_grid = lambda_grid * res_var_scale,
                                                          lambda_ratio = lambda_ratio)
          lambda_signal <- lambda_ratio * lambda_parallel
        } else {
          warning("gcv_lambda = TRUE but projection_mats$P_Noise is missing")
        }
      }

      n_reg <- ncol(X_eff)
      K_penalty_mat <- matrix(0, n_reg, n_reg)
      if (!is.null(projection_mats$P_GD)) {
        K_penalty_mat <- K_penalty_mat + lambda_gd * projection_mats$P_GD
      }
      if (!is.null(projection_mats$P_Unique)) {
        K_penalty_mat <- K_penalty_mat + lambda_unique * projection_mats$P_Unique
      }
      if (!is.null(projection_mats$P_Noise)) {
        K_penalty_mat <- K_penalty_mat + lambda_parallel * projection_mats$P_Noise
      }
      if (!is.null(projection_mats$P_Signal)) {
        K_penalty_mat <- K_penalty_mat + lambda_signal * projection_mats$P_Signal
      }
    }
    list(K_penalty_diag = NULL, K_penalty_mat = K_penalty_mat)
  } else {
    if (is.null(K_penalty_diag)) {
      lambda_parallel <- lambda_values$lambda_parallel %||% 1.0
      lambda_signal   <- lambda_values$lambda_perp_signal %||% (0.05 * lambda_parallel)
      lambda_gd       <- lambda_values$lambda_gd %||% lambda_parallel
      lambda_unique   <- lambda_values$lambda_unique %||% lambda_parallel

      if (gcv_lambda) {
        if (!is.null(projection_mats$P_Noise)) {
          lambda_ratio <- lambda_signal / lambda_parallel
          lambda_parallel <- ndx_gcv_tune_lambda_parallel(Y_eff, X_eff,
                                                          projection_mats$P_Noise,
                                                          lambda_grid = lambda_grid * res_var_scale,
                                                          lambda_ratio = lambda_ratio)
          lambda_signal <- lambda_ratio * lambda_parallel
        } else {
          warning("gcv_lambda = TRUE but projection_mats$P_Noise is missing")
        }
      }

      n_reg <- ncol(X_eff)
      K_penalty_diag <- rep(0, n_reg)
      if (!is.null(projection_mats$P_GD)) {
        K_penalty_diag <- K_penalty_diag + lambda_gd * diag(projection_mats$P_GD)
      }
      if (!is.null(projection_mats$P_Unique)) {
        K_penalty_diag <- K_penalty_diag + lambda_unique * diag(projection_mats$P_Unique)
      }
      if (!is.null(projection_mats$P_Noise)) {
        K_penalty_diag <- K_penalty_diag + lambda_parallel * diag(projection_mats$P_Noise)
      }
      if (!is.null(projection_mats$P_Signal)) {
        K_penalty_diag <- K_penalty_diag + lambda_signal * diag(projection_mats$P_Signal)
      }
    }
    list(K_penalty_diag = K_penalty_diag, K_penalty_mat = NULL)
  }
}

#' Solve Weighted Ridge System
#'
#' Internal helper performing the Cholesky solve for anisotropic ridge.
#'
#' @keywords internal
#' @return List with elements `betas` and `chol_RtR` (NULL when weights are a matrix).
solve_weighted_ridge <- function(Y_eff, X_eff, K_penalty_diag = NULL,
                                 K_penalty_mat = NULL, W_eff = NULL,
                                 use_penalty_matrix = FALSE,
                                 x_names = NULL, y_names = NULL) {

  col_vars <- matrixStats::colVars(X_eff, na.rm = TRUE)
  diag_pen <- if (use_penalty_matrix) diag(K_penalty_mat) else K_penalty_diag
  if (any(col_vars < .Machine$double.eps^0.5 &
          diag_pen[which(col_vars < .Machine$double.eps^0.5)] < .Machine$double.eps^0.5)) {
    warning("One or more regressors in X_eff have near-zero variance AND near-zero penalty. Solution may be unstable.")
  } else if (any(col_vars < .Machine$double.eps^0.5)) {
    warning("One or more regressors in X_eff have near-zero variance (but may be stabilized by penalty).")
  }

  if (is.null(W_eff)) {
    XtX <- crossprod(X_eff)
    XtY <- crossprod(X_eff, Y_eff)
  } else if (is.matrix(W_eff)) {
    betas <- matrix(NA_real_, ncol = ncol(Y_eff), nrow = ncol(X_eff))
    for (v in seq_len(ncol(Y_eff))) {
      w_sqrt <- sqrt(W_eff[, v])
      Xw <- X_eff * w_sqrt
      Yw <- Y_eff[, v] * w_sqrt
      XtX <- crossprod(Xw)
      XtY <- crossprod(Xw, Yw)

      if (use_penalty_matrix) {
        K_mat_safe <- K_penalty_mat
        diag(K_mat_safe) <- pmax(diag(K_mat_safe), .Machine$double.eps)
        lhs <- XtX + K_mat_safe
      } else {
        safe_K_penalty_diag <- pmax(K_penalty_diag, .Machine$double.eps)
        lhs <- XtX
        diag(lhs) <- diag(lhs) + safe_K_penalty_diag
      }

      b <- tryCatch({
        chol_decomp <- chol(lhs)
        chol2inv(chol_decomp) %*% XtY
      }, error = function(e) NULL)
      if (!is.null(b)) betas[, v] <- b
    }
    if (!is.null(x_names)) rownames(betas) <- x_names
    if (!is.null(y_names)) colnames(betas) <- y_names
    return(list(betas = betas, chol_RtR = NULL))
  } else {
    w_sqrt <- sqrt(W_eff)
    Xw <- X_eff * w_sqrt
    Yw <- Y_eff * w_sqrt
    XtX <- crossprod(Xw)
    XtY <- crossprod(Xw, Yw)
  }

  if (ncol(X_eff) == 0) {
    warning("X_eff has zero columns after NA removal. Cannot solve.")
    return(NULL)
  }

  if (use_penalty_matrix) {
    K_mat_safe <- K_penalty_mat
    diag(K_mat_safe) <- pmax(diag(K_mat_safe), .Machine$double.eps)
    lhs <- XtX + K_mat_safe
  } else {
    safe_K_penalty_diag <- pmax(K_penalty_diag, .Machine$double.eps)
    lhs <- XtX
    diag(lhs) <- diag(lhs) + safe_K_penalty_diag
  }

  betas <- NULL
  chol_decomp <- NULL
  tryCatch({
    chol_decomp <- chol(lhs)
    betas <- chol2inv(chol_decomp) %*% XtY
  }, error = function(e) {
    warning(paste("Solving anisotropic ridge regression failed (Cholesky method):", e$message))
    betas <<- NULL
    chol_decomp <<- NULL
  })

  if (is.null(betas)) {
    return(NULL)
  }

  if (!is.null(x_names)) rownames(betas) <- x_names
  if (!is.null(y_names)) colnames(betas) <- y_names

  list(betas = betas, chol_RtR = chol_decomp)
}



#' Solve Anisotropic Ridge Regression Problem
#'
#' Solves for beta coefficients in a ridge regression model: Y = X beta + E,
#' minimizing ||Y - X beta||^2 + beta^T K beta, where K is a diagonal penalty matrix.
#' Handles missing data specified by `na_mask` by removing corresponding timepoints.
#' This solver underpins the modular anisotropic ridge API.
#'
#' @param Y_whitened A numeric matrix of whitened dependent variables (timepoints x voxels/responses).
#' @param X_whitened A numeric matrix of whitened regressors (timepoints x n_regressors).
#' @param K_penalty_diag A numeric vector containing the diagonal elements of the penalty matrix K.
#'   Its length must be equal to `ncol(X_whitened)`.
#' @param na_mask Optional. A logical vector where TRUE indicates timepoints to
#'   exclude. If NULL, all timepoints are used.
#' @param projection_mats Optional list of projection matrices produced by
#'   `ndx_compute_projection_matrices`. If provided along with `lambda_values`,
#'   `K_penalty_diag` (or `K_penalty_mat` when `use_penalty_matrix = TRUE`) is
#'   constructed automatically.
#' @param lambda_values Optional list with elements `lambda_parallel`,
#'   `lambda_perp_signal`, `lambda_gd`, and `lambda_unique` used when a penalty
#'   matrix or diagonal is not given.
#' @param K_penalty_mat Optional full penalty matrix to use instead of
#'   `K_penalty_diag`. When supplied, it must be a square matrix with dimension
#'   equal to `ncol(X_whitened)`.
#' @param use_penalty_matrix Logical, when TRUE a full penalty matrix is used.
#'   The matrix is either supplied via `K_penalty_mat` or constructed from
#'   `projection_mats` and `lambda_values`. The default FALSE preserves the
#'   historic diagonal-only behaviour.
#' @param weights Optional numeric vector or matrix of weights to apply to each
#'   timepoint. If a matrix, its dimensions must match `Y_whitened` and the
#'   weighting is applied per voxel.
#' @param gcv_lambda Logical, if TRUE perform GCV tuning of `lambda_parallel`
#'   (and by extension `lambda_perp_signal`) using `ndx_gcv_tune_lambda_parallel`.
#' @param res_var_scale Numeric, residual variance estimate used to scale the
#'   lambda search grid when `gcv_lambda` is TRUE.
#' @param lambda_grid Numeric vector of candidate lambda values for GCV when
#'   `gcv_lambda` is TRUE. Defaults to `10^seq(-2, 2, length.out = 5)`.
#' @return A list with elements:
#'   \describe{
#'     \item{betas}{Matrix of estimated coefficients (n_regressors x voxels).}
#'     \item{chol_RtR}{Cholesky factor of the penalized normal equations. NULL when weights are a matrix.}
#'   }
#'   Returns NULL if inputs are invalid or the problem cannot be solved.
#' @export ndx_solve_anisotropic_ridge
ndx_solve_anisotropic_ridge <- function(Y_whitened, X_whitened, K_penalty_diag = NULL,
                                        na_mask = NULL, projection_mats = NULL,
                                        lambda_values = NULL, weights = NULL,
                                        gcv_lambda = FALSE, res_var_scale = 1.0,
                                        lambda_grid = 10^seq(-2, 2, length.out = 5),
                                        K_penalty_mat = NULL,
                                        use_penalty_matrix = FALSE) {

  if (!validate_ridge_inputs(Y_whitened, X_whitened, K_penalty_diag, na_mask,
                             projection_mats, lambda_values, weights,
                             K_penalty_mat, use_penalty_matrix)) {
    return(NULL)
  }

  na_res <- apply_na_mask_and_weights(Y_whitened, X_whitened, na_mask, weights)
  if (is.null(na_res)) return(NULL)
  Y_eff <- na_res$Y
  X_eff <- na_res$X
  W_eff <- na_res$W

  pen <- construct_penalty(X_eff, Y_eff, K_penalty_diag,
                           projection_mats, lambda_values,
                           gcv_lambda, res_var_scale, lambda_grid,
                           K_penalty_mat, use_penalty_matrix)
  K_penalty_diag <- pen$K_penalty_diag
  K_penalty_mat  <- pen$K_penalty_mat

  solve_weighted_ridge(Y_eff, X_eff,
                       K_penalty_diag = K_penalty_diag,
                       K_penalty_mat = K_penalty_mat,
                       W_eff = W_eff,
                       use_penalty_matrix = use_penalty_matrix,
                       x_names = colnames(X_whitened),
                       y_names = colnames(Y_whitened))
}

#' Extract Task-Related Beta Coefficients
#'
#' Extracts beta coefficients corresponding to specified task regressors from the full
#' set of beta coefficients. For Sprint 1, this function primarily handles extraction.
#' The "unwhitening" aspect is a placeholder for future refinement.
#' This helper is part of the modular anisotropic ridge API.
#'
#' @param betas_whitened A matrix (total_regressors x voxels) of whitened beta coefficients,
#'   typically the output of `ndx_solve_ridge`.
#' @param X_whitened_colnames A character vector of column names for the whitened design matrix (`X_whitened`)
#'   that was used to estimate `betas_whitened`. This is used to identify task regressors by name.
#' @param task_regressor_names A character vector containing the names of the task regressors
#'   to be extracted.
#' @param ar_coeffs_global Optional. The global AR coefficients (vector of length `order`) that were used
#'   to whiten the design matrix `X`. If provided, a conceptual note about unwhitening can be considered.
#'   (Currently not used for actual unwhitening).
#' @param unwhiten Logical, if TRUE, attempt to unwhiten betas. Default is FALSE (not yet implemented).
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
ndx_extract_task_betas <- function(betas_whitened, X_whitened_colnames, task_regressor_names, 
                                   ar_coeffs_global = NULL, unwhiten = FALSE) {

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
  
  if (unwhiten) {
    warning("Unwhitening of betas (argument unwhiten=TRUE) is not yet implemented. Returning whitened betas.")
    # Placeholder for future unwhitening logic using ar_coeffs_global if provided
  }
  
  return(extracted_betas)
} 