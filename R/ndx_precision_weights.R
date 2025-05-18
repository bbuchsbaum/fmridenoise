#' Compute Precision Weights from RPCA Sparse Matrix
#'
#' Given the sparse component matrix `S_matrix` from RPCA, this helper
#' constructs precision weights that down-weight timepoints with large
#' spike amplitudes.
#'
#' @param S_matrix Numeric matrix of RPCA sparse components (timepoints x voxels).
#' @param mad_floor Numeric minimum value for the MAD scaling. Default 1e-6.
#' @return Numeric matrix of precision weights with the same dimensions as
#'   `S_matrix`.
#' @export ndx_precision_weights_from_S
ndx_precision_weights_from_S <- function(S_matrix, mad_floor = 1e-6) {
  if (!is.matrix(S_matrix) || !is.numeric(S_matrix)) {
    stop("S_matrix must be a numeric matrix.")
  }
  if (!is.numeric(mad_floor) || length(mad_floor) != 1 || mad_floor <= 0) {
    stop("mad_floor must be a positive numeric scalar.")
  }

  mad_S <- max(stats::mad(abs(S_matrix), na.rm = TRUE), mad_floor)
  1 / (abs(S_matrix) / mad_S + 1)^2
}
