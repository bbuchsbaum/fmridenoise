#' @keywords internal
# NO @export for prepare_hrf_response_data
prepare_hrf_response_data <- function(Y_fmri, pass0_residuals, run_idx,
                                      validated_spike_TR_mask, user_options) {
  R2_pass0_vox <- calculate_R2_voxelwise(Y_fmri, pass0_residuals)
  good_voxels_idx <- which(R2_pass0_vox >= user_options$good_voxel_R2_threshold)
  n_good_voxels <- length(good_voxels_idx)

  if (n_good_voxels < user_options$hrf_min_good_voxels) {
    warning(sprintf(
      "Insufficient good voxels found (%d) based on R2 threshold (%.2f). Minimum required: %d. Skipping HRF estimation.",
      n_good_voxels, user_options$good_voxel_R2_threshold,
      user_options$hrf_min_good_voxels
    ))
    return(NULL)
  }

  Y_good_voxels <- Y_fmri[, good_voxels_idx, drop = FALSE]
  if (ncol(Y_good_voxels) == 0) {
    warning(
      "No good voxels selected (ncol=0) despite passing minimum count. Check R2 threshold or data. Skipping HRF estimation."
    )
    return(NULL)
  }

  valid_TRs_for_hrf_estimation_mask <- !validated_spike_TR_mask
  if (sum(valid_TRs_for_hrf_estimation_mask) == 0) {
    warning(
      "No valid TRs available for HRF estimation after applying spike_TR_mask. Skipping HRF estimation."
    )
    return(NULL)
  }

  ybar_all_trs <- matrixStats::rowMedians(Y_good_voxels, na.rm = TRUE)
  ybar_clean <- ybar_all_trs[valid_TRs_for_hrf_estimation_mask]
  block_ids_for_cv <- run_idx[valid_TRs_for_hrf_estimation_mask]

  return(list(
    ybar_clean = ybar_clean,
    block_ids_for_cv = block_ids_for_cv,
    valid_TRs_for_hrf_estimation_mask = valid_TRs_for_hrf_estimation_mask,
    good_voxels_indices = good_voxels_idx,
    Y_good_voxels_all_trs = Y_good_voxels,
    n_good_voxels = n_good_voxels
  ))
}
