context("HRF clustering utilities")

# Test .perform_hrf_clustering

test_that(".perform_hrf_clustering returns expected clusters", {
  set.seed(123)
  n_time <- 30
  t <- seq(0, 2*pi, length.out = n_time)
  pattern1 <- sin(t)
  pattern2 <- cos(t)

  Y <- cbind(
    pattern1 + matrix(rnorm(n_time*3, sd = 0.01), n_time, 3),
    pattern2 + matrix(rnorm(n_time*3, sd = 0.01), n_time, 3)
  )

  user_opts <- list(hrf_min_good_voxels = 1L)

  cl_res <- ndx:::.perform_hrf_clustering(
    Y_for_clustering = Y,
    num_clusters = 2L,
    user_options = user_opts,
    verbose = FALSE
  )

  expect_equal(cl_res$num_clusters_eff, 2L)
  expect_equal(length(cl_res$voxel_cluster_ids), ncol(Y))
  expect_length(unique(cl_res$voxel_cluster_ids[1:3]), 1L)
  expect_length(unique(cl_res$voxel_cluster_ids[4:6]), 1L)
  expect_false(cl_res$voxel_cluster_ids[1] == cl_res$voxel_cluster_ids[4])
})

# Test .auto_adapt_hrf_clusters

test_that(".auto_adapt_hrf_clusters merges and reassigns clusters", {
  set.seed(456)
  n_time <- 30
  t <- seq(0, 2*pi, length.out = n_time)
  pattern_a <- sin(t)
  pattern_b <- sin(t) + rnorm(n_time, sd = 0.001)  # Highly correlated with pattern_a
  pattern_c <- cos(t)

  Y <- cbind(
    pattern_a + matrix(rnorm(n_time*4, sd = 0.01), n_time, 4),
    pattern_b + matrix(rnorm(n_time*4, sd = 0.01), n_time, 4),
    pattern_c + matrix(rnorm(n_time*2, sd = 0.01), n_time, 2)
  )

  initial_ids <- c(rep(1L,4), rep(2L,4), rep(3L,2))

  adapt_res <- ndx:::.auto_adapt_hrf_clusters(
    Y_for_clustering = Y,
    cluster_ids = initial_ids,
    merge_corr_thresh = 0.95,
    min_cluster_size = 3L,
    verbose = FALSE
  )

  expect_equal(adapt_res$num_clusters_eff, 1L)
  expect_equal(length(adapt_res$voxel_cluster_ids), ncol(Y))
  expect_equal(unique(adapt_res$voxel_cluster_ids), 1L)
})
