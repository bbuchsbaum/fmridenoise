context("ndx_run_gdlite")

test_that("run_gdlite returns expected structure on toy data", {
  set.seed(1)
  n_time <- 6
  run_idx <- rep(1:2, each = 3)
  TR <- 1.0
  Y <- matrix(rnorm(n_time * 4), n_time, 4)

  events <- data.frame(
    onsets = c(0, 3),
    durations = c(1, 1),
    condition = "A",
    blockids = 1:2
  )

  res <- ndx_run_gdlite(
    Y_fmri = Y,
    events = events,
    run_idx = run_idx,
    TR = TR,
    poly_degree = 0,
    k_max = 2,
    r2_thresh_noise_pool = 0,
    tsnr_thresh_noise_pool = -Inf,
    r2_thresh_good_voxels = 0,
    detrend_for_tsnr = FALSE,
    perform_final_glm = FALSE
  )

  expect_true(is.list(res))
  expect_true(all(c("selected_pcs", "K_star", "X_gd") %in% names(res)))
  if (res$K_star > 0) {
    expect_true(is.matrix(res$selected_pcs))
    expect_equal(nrow(res$selected_pcs), n_time)
    expect_equal(ncol(res$selected_pcs), res$K_star)
  } else {
    expect_null(res$selected_pcs)
  }
  expect_true(is.matrix(res$X_gd))
})
