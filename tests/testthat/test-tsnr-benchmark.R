context("Benchmark calculate_tsnr vectorization")

# Old implementation with per-voxel loop
calculate_tsnr_loop <- function(Y_data, detrend = FALSE, robust = FALSE) {
  n_voxels <- ncol(Y_data)
  n_timepoints <- nrow(Y_data)
  tsnr_values <- numeric(n_voxels)
  time_vector <- seq_len(n_timepoints)
  for (v_idx in seq_len(n_voxels)) {
    voxel_series <- Y_data[, v_idx]
    valid_mask <- !is.na(voxel_series)
    if (sum(valid_mask) < 2) {
      tsnr_values[v_idx] <- NA_real_
      next
    }
    voxel_series_valid <- voxel_series[valid_mask]
    time_vector_valid <- time_vector[valid_mask]
    if (detrend) {
      if (length(unique(time_vector_valid)) > 1) {
        fit <- tryCatch(stats::lm(voxel_series_valid ~ time_vector_valid), error = function(e) NULL)
        if (!is.null(fit)) {
          voxel_series_processed <- stats::residuals(fit)
        } else {
          voxel_series_processed <- voxel_series_valid - mean(voxel_series_valid, na.rm = TRUE)
        }
      } else {
        voxel_series_processed <- voxel_series_valid - mean(voxel_series_valid, na.rm = TRUE)
      }
    } else {
      voxel_series_processed <- voxel_series_valid
    }
    if (robust) {
      mean_val <- stats::median(voxel_series_processed, na.rm = TRUE)
      sd_val <- stats::mad(voxel_series_processed, na.rm = TRUE)
    } else {
      mean_val <- mean(voxel_series_processed, na.rm = TRUE)
      sd_val <- stats::sd(voxel_series_processed, na.rm = TRUE)
    }
    if (is.na(sd_val) || sd_val < 1e-9) {
      tsnr_values[v_idx] <- NA_real_
    } else {
      tsnr_values[v_idx] <- mean_val / sd_val
    }
  }
  tsnr_values
}

test_that("vectorized calculate_tsnr is faster than loop version", {
  testthat::skip_if_not_installed("microbenchmark")
  set.seed(123)
  Y <- matrix(rnorm(200 * 100), nrow = 200, ncol = 100)
  bench <- microbenchmark::microbenchmark(
    loop = calculate_tsnr_loop(Y),
    vectorized = calculate_tsnr(Y),
    times = 5L
  )
  expect_lt(median(bench$time[bench$expr == "vectorized"]),
            median(bench$time[bench$expr == "loop"]))
})
