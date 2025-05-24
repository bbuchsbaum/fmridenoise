#' Compute Precision Weights from RPCA Sparse Matrix
#'
#' Given the sparse component matrix `S_matrix` from RPCA, this helper
#' constructs precision weights that down-weight timepoints with large
#' spike amplitudes.
#'
#' @param S_matrix Numeric matrix of RPCA sparse components (timepoints x voxels).
#' @param mad_floor Numeric minimum value for the MAD scaling. Default 1e-6.
#' @param na_zero Logical; if `TRUE`, any `NA` or infinite values in
#'   `S_matrix` are treated as zero (yielding a weight of one). Default `TRUE`.
#' @return Numeric matrix of precision weights with the same dimensions as
#'   `S_matrix`.
#' @export ndx_precision_weights_from_S
ndx_precision_weights_from_S <- function(S_matrix, mad_floor = 1e-6,
                                         na_zero = TRUE) {
  if (!is.matrix(S_matrix) || !is.numeric(S_matrix)) {
    stop("S_matrix must be a numeric matrix.")
  }
  if (!is.numeric(mad_floor) || length(mad_floor) != 1 || mad_floor <= 0) {
    stop("mad_floor must be a positive numeric scalar.")
  }
  if (!is.logical(na_zero) || length(na_zero) != 1) {
    stop("na_zero must be a single logical value.")
  }

  abs_S <- abs(S_matrix)
  if (na_zero) {
    abs_S[!is.finite(abs_S)] <- 0
  }

  mad_S <- max(stats::mad(abs_S, na.rm = TRUE), mad_floor)
  weights <- 1 / (abs_S / mad_S + 1)^2
  if (na_zero) {
    weights[!is.finite(weights)] <- 1
  }
  weights
}
