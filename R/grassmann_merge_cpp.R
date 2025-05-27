#' Rcpp implementation of iterative Grassmann merge
#'
#' @param V_list_valid List of matrices (voxels x k_r).
#' @param k_target_global Target number of global components.
#' @return Matrix V_global with attribute `svd_times`.
#' @keywords internal
.grassmann_merge_iterative_rcpp <- function(V_list_valid, k_target_global) {
  res <- .grassmann_merge_iterative_cpp(V_list_valid, as.integer(k_target_global))
  Vg <- res$V_global
  attr(Vg, "svd_times") <- res$svd_times
  Vg
}
