#' Perform K-Medoids Clustering on Voxel Time Courses
#' @keywords internal
.perform_hrf_clustering <- function(Y_for_clustering, num_clusters, user_options, verbose) {
  # Y_for_clustering is timepoints x good_voxels_for_clustering
  # We cluster voxels, so pam needs voxels x timepoints, or dist on voxels
  if (verbose) message(sprintf("  Performing K-Medoids clustering for %d clusters on %d voxels.", num_clusters, ncol(Y_for_clustering)))
  
  if (num_clusters == 1 || ncol(Y_for_clustering) <= num_clusters || ncol(Y_for_clustering) < (user_options$hrf_min_good_voxels %||% 10) ) {
    if (verbose && num_clusters > 1) message("    Number of good voxels too small for requested clusters, or num_clusters=1. Assigning all to cluster 1.")
    return(list(voxel_cluster_ids = rep(1L, ncol(Y_for_clustering)), 
                num_clusters_eff = 1L, 
                medoid_indices = if(ncol(Y_for_clustering)>0) 1L else integer(0) )) # Default to first voxel as medoid if any
  }
  
  # cluster::pam expects samples in rows, features in columns. 
  # To cluster voxels (features here), we can transpose Y_for_clustering.
  # Or, supply a distance matrix. Transposing is simpler if data fits memory.
  # Consider if scaling/normalizing voxels before clustering is needed (e.g. z-score each voxel's timecourse)
  pam_fit <- NULL
  tryCatch({
    # Note: pam can be slow for large N (voxels) and large P (timepoints per voxel)
    # For very large number of voxels, might need to subsample or use clara.
    # Data for pam should be n_observations (voxels) x n_variables (timepoints)
    pam_fit <- cluster::pam(t(Y_for_clustering), k = num_clusters, diss = FALSE, metric="euclidean", stand = FALSE)
  }, error = function(e) {
    warning(sprintf("K-Medoids clustering (pam) failed: %s. Assigning all to cluster 1.", e$message))
    pam_fit <- NULL # Keep pam_fit local on error
  })
  
  if (is.null(pam_fit)) {
    return(list(voxel_cluster_ids = rep(1L, ncol(Y_for_clustering)), 
                num_clusters_eff = 1L, 
                medoid_indices = if(ncol(Y_for_clustering)>0) 1L else integer(0) ))
  }
  
  return(list(voxel_cluster_ids = pam_fit$clustering,
              num_clusters_eff = length(unique(pam_fit$clustering)),
              medoid_indices = pam_fit$id.med # Indices of medoids *within the input to pam* (i.e., among good_voxels_idx)
              ))
}

#' Auto-adapt HRF clusters by merging highly similar or tiny clusters
#'
#' This simplified routine merges clusters whose mean time courses are
#' highly correlated and reassigns clusters with very few voxels to the
#' most similar larger cluster.
#' @keywords internal
.auto_adapt_hrf_clusters <- function(Y_for_clustering, cluster_ids,
                                    merge_corr_thresh = 0.95,
                                    min_cluster_size = 5L,
                                    verbose = FALSE) {
  if (is.null(cluster_ids) || length(unique(cluster_ids)) <= 1) {
    return(list(voxel_cluster_ids = cluster_ids,
                num_clusters_eff = length(unique(cluster_ids)),
                medoid_indices = if(length(cluster_ids) > 0) 1L else integer(0)))
  }

  cluster_ids <- as.integer(cluster_ids)
  changed <- TRUE
  while (changed) {
    changed <- FALSE
    cl_levels <- sort(unique(cluster_ids))
    if (length(cl_levels) <= 1) break

    cluster_means <- sapply(cl_levels, function(cl)
      rowMeans(Y_for_clustering[, cluster_ids == cl, drop = FALSE]))
    corr_mat <- stats::cor(cluster_means)
    diag(corr_mat) <- NA
    high_pairs <- which(corr_mat > merge_corr_thresh, arr.ind = TRUE)
    if (nrow(high_pairs) > 0) {
      pair <- high_pairs[1, ]
      cl1 <- cl_levels[pair[1]]; cl2 <- cl_levels[pair[2]]
      sz1 <- sum(cluster_ids == cl1); sz2 <- sum(cluster_ids == cl2)
      if (verbose)
        message(sprintf("  Merging clusters %d and %d due to high correlation (%.2f)",
                        cl1, cl2, corr_mat[pair[1], pair[2]]))
      keep_cl <- if (sz1 >= sz2) cl1 else cl2
      drop_cl <- if (sz1 >= sz2) cl2 else cl1
      cluster_ids[cluster_ids == drop_cl] <- keep_cl
      changed <- TRUE
    }
  }

  # Reassign small clusters
  cl_levels <- sort(unique(cluster_ids))
  if (length(cl_levels) > 1) {
    cluster_means <- sapply(cl_levels, function(cl)
      rowMeans(Y_for_clustering[, cluster_ids == cl, drop = FALSE]))
    for (cl in cl_levels) {
      vox_idx <- which(cluster_ids == cl)
      if (length(vox_idx) < min_cluster_size) {
        other_cls <- setdiff(cl_levels, cl)
        if (length(other_cls) == 0) next
        cor_to_others <- stats::cor(rowMeans(Y_for_clustering[, vox_idx, drop = FALSE]),
                                    cluster_means[, match(other_cls, cl_levels), drop = FALSE])
        nearest <- other_cls[which.max(cor_to_others)]
        if (verbose)
          message(sprintf("  Reassigning small cluster %d (n=%d) to %d", cl, length(vox_idx), nearest))
        cluster_ids[vox_idx] <- nearest
      }
    }
  }

  # Renumber sequentially
  final_levels <- sort(unique(cluster_ids))
  mapping <- setNames(seq_along(final_levels), final_levels)
  cluster_ids <- mapping[as.character(cluster_ids)]
  final_levels <- sort(unique(cluster_ids))

  medoid_indices <- integer(length(final_levels))
  for (i in seq_along(final_levels)) {
    cl <- final_levels[i]
    vox_idx <- which(cluster_ids == cl)
    if (length(vox_idx) == 1) {
      medoid_indices[i] <- vox_idx
    } else {
      pam_res <- cluster::pam(t(Y_for_clustering[, vox_idx, drop = FALSE]), k = 1)
      medoid_indices[i] <- vox_idx[pam_res$id.med]
    }
  }

  list(voxel_cluster_ids = cluster_ids,
       num_clusters_eff = length(final_levels),
       medoid_indices = medoid_indices)
}
