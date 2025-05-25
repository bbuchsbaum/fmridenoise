#' Helper: Robust PCA on a Single Run
#'
#' Performs RPCA on a single run's residual matrix and extracts voxel-space
#' components along with diagnostic information.
#'
#' @param E_run Numeric matrix of residuals (timepoints x voxels).
#' @param k_target Integer target rank for the low-rank component.
#' @param opts List of options. Recognized elements: `rpca_lambda_auto`,
#'   `rpca_lambda_fixed`, `rpca_mu`, `rpca_term_delta`, `rpca_max_iter`,
#'   `rpca_trace`, `spike_mad_thresh`, and `spike_percentile_thresh`.
#' @return List with elements `V_r`, `glitch_ratio`, `S_matrix_TxV`,
#'   and `spike_TR_mask`. Returns basic empty structures if input is invalid.
#' @keywords internal
.ndx_rpca_single_run <- function(E_run, k_target = 1L, opts = list()) {
  if (!is.matrix(E_run) || nrow(E_run) == 0 || ncol(E_run) == 0) {
    return(list(
      V_r = NULL,
      glitch_ratio = NA_real_,
      S_matrix_TxV = matrix(numeric(0), nrow(E_run), ncol(E_run)),
      spike_TR_mask = logical(nrow(E_run))
    ))
  }

  defaults <- list(
    rpca_lambda_auto = TRUE,
    rpca_lambda_fixed = NULL,
    rpca_mu = NULL,
    rpca_term_delta = 1e-6,
    rpca_max_iter = 2000,
    rpca_trace = FALSE,
    spike_mad_thresh = 3.0,
    spike_percentile_thresh = 0.98
  )
  opts <- utils::modifyList(defaults, opts)

  E_t <- t(E_run)
  lambda_r <- if (opts$rpca_lambda_auto) {
    1 / sqrt(max(dim(E_t)))
  } else {
    if (is.null(opts$rpca_lambda_fixed))
      stop("rpca_lambda_fixed must be provided when rpca_lambda_auto is FALSE")
    opts$rpca_lambda_fixed
  }

  rpca_args <- list(
    M = E_t,
    lambda = lambda_r,
    term.delta = opts$rpca_term_delta,
    max.iter = opts$rpca_max_iter,
    trace = opts$rpca_trace
  )
  if (!is.null(opts$rpca_mu)) rpca_args$mu <- opts$rpca_mu

  rp <- tryCatch({
    do.call(rpca::rpca, rpca_args)
  }, error = function(e) NULL)

  if (is.null(rp) || is.null(rp$L)) {
    return(list(
      V_r = NULL,
      glitch_ratio = NA_real_,
      S_matrix_TxV = matrix(0, nrow(E_run), ncol(E_run)),
      spike_TR_mask = rep(FALSE, nrow(E_run))
    ))
  }

  L_t <- rp$L
  k_eff <- min(k_target, nrow(L_t), ncol(L_t))
  svd_L <- svd(L_t, nu = k_eff, nv = 0)
  V_r <- svd_L$u

  S_t <- rp$S
  e_L <- sum(L_t^2)
  e_S <- sum(S_t^2)
  glitch_ratio <- if (e_L > 0) e_S / e_L else NA_real_

  if (!is.null(S_t) && length(S_t) > 0) {
    s_star <- apply(abs(S_t), 2, stats::median)
    mad_s <- stats::mad(s_star, constant = 1)
    if (!is.finite(mad_s) || mad_s < 1e-9) {
      thr <- stats::quantile(s_star, probs = opts$spike_percentile_thresh)
    } else {
      thr <- median(s_star) + opts$spike_mad_thresh * mad_s
    }
    spike_mask <- s_star > thr
  } else {
    spike_mask <- rep(FALSE, nrow(E_run))
  }

  list(
    V_r = V_r,
    glitch_ratio = glitch_ratio,
    S_matrix_TxV = t(S_t),
    spike_TR_mask = spike_mask
  )
}

#' Helper: Merge Voxel Components
#'
#' Combines a list of voxel-space component matrices into a global basis.
#'
#' @param V_list List of matrices (voxels x k).
#' @param k_target Desired number of global components.
#' @param strategy "concat_svd" or "iterative".
#' @return Matrix `V_global` or `NULL` if merging fails.
#' @keywords internal
.ndx_merge_voxel_components <- function(V_list, k_target, 
                                        strategy = c("concat_svd", "iterative")) {
  strategy <- match.arg(strategy)
  V_list_valid <- Filter(function(x) is.matrix(x) && ncol(x) > 0, V_list)
  if (length(V_list_valid) == 0 || k_target <= 0) return(NULL)

  if (strategy == "concat_svd") {
    V_all <- do.call(cbind, V_list_valid)
    if (is.null(V_all) || ncol(V_all) == 0) return(NULL)
    k_for_svd <- min(k_target, nrow(V_all), ncol(V_all))
    svd_res <- svd(V_all, nu = k_for_svd, nv = 0)
    return(svd_res$u)
  }

  .grassmann_merge_iterative(V_list_valid, k_target)
}


