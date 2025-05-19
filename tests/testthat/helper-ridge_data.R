# Helper functions for ridge regression tests

# Generates simulated ridge regression data with optional NA periods.
# Returns a list containing Y, X, true_betas, and na_mask.
.create_ridge_test_data <- function(n_timepoints, n_regressors, n_voxels,
                                    true_betas, noise_sd = 1.0, na_head = 0) {
  set.seed(123)  # Ensure reproducibility across tests
  X <- matrix(rnorm(n_timepoints * n_regressors), n_timepoints, n_regressors)

  if (n_regressors > 0) X[, 1] <- 1
  if (n_regressors > 1) X[, 2] <- seq(0, 1, length.out = n_timepoints)
  if (n_regressors > 2) X[, 3] <- (seq(0, 1, length.out = n_timepoints))^2
  if (n_regressors > 3) X[, 4] <- sin(seq(0, 2 * pi, length.out = n_timepoints))

  if (n_regressors > 0) {
    if (is.vector(true_betas)) {
      if (length(true_betas) != n_regressors)
        stop("Length of true_betas vector must match n_regressors")
      true_betas_matrix <- matrix(rep(true_betas, n_voxels),
                                  nrow = n_regressors, ncol = n_voxels)
    } else if (is.matrix(true_betas)) {
      if (nrow(true_betas) != n_regressors || ncol(true_betas) != n_voxels)
        stop("Dimensions of true_betas matrix are incorrect.")
      true_betas_matrix <- true_betas
    } else {
      stop("true_betas must be a vector or a matrix.")
    }
    Y_signal <- X %*% true_betas_matrix
  } else {
    true_betas_matrix <- matrix(0, nrow = 0, ncol = n_voxels)
    Y_signal <- matrix(0, nrow = n_timepoints, ncol = n_voxels)
  }

  Y_noise <- matrix(rnorm(n_timepoints * n_voxels, sd = noise_sd),
                    n_timepoints, n_voxels)
  Y <- Y_signal + Y_noise

  na_mask_vec <- rep(FALSE, n_timepoints)
  if (na_head > 0 && n_timepoints > na_head) {
    Y[1:na_head, ] <- NA
    na_mask_vec[1:na_head] <- TRUE
  }

  if (n_regressors > 0) colnames(X) <- paste0("Reg", 1:n_regressors)
  if (n_voxels > 0) colnames(Y) <- paste0("Voxel", 1:n_voxels)

  list(Y = Y, X = X, true_betas = true_betas_matrix, na_mask = na_mask_vec)
}

