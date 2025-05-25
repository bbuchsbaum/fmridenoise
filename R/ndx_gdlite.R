#' @title Extract GLMdenoise-Lite Principal Components
#'
#' @description
#' Performs a basic GLM to obtain residuals from fMRI data after regressing out
#' a base set of nuisance regressors (e.g., task, motion, trends). Then, PCA is
#' applied to these residuals to extract temporal principal components, similar to
#' a simplified GLMdenoise approach.
#'
#' @param Y_fmri Numeric matrix of fMRI data (timepoints x voxels).
#' @param X_nuisance_base Numeric matrix representing the design matrix of
#'   effects to regress out before PCA (e.g., task effects from a canonical HRF,
#'   motion parameters, polynomial trends). Must have the same number of rows
#'   as `Y_fmri`.
#' @param n_pcs Integer, the number of principal components to extract.
#' @param voxel_mask Optional. A logical or numeric vector indicating which voxels
#'   (columns of `Y_fmri`) to include in the PCA. If logical, length must match
#'   `ncol(Y_fmri)`. If numeric, contains column indices. If `NULL` (default),
#'   PCA will be performed on a subset of high-variance voxels (see `min_variance_voxels`, `max_variance_voxels_pca_prop`).
#' @param min_variance_voxels Integer, the minimum number of voxels to use for PCA
#'   if `voxel_mask` is `NULL`. Defaults to `max(100, n_pcs * 5)`.
#' @param max_variance_voxels_pca_prop Numeric, proportion of total voxels to consider as
#'   candidates for high-variance selection if `voxel_mask` is `NULL`. E.g. 0.2 means top 20% variance voxels.
#'   The actual number used will be `min(max_variance_voxels, number_of_voxels_meeting_prop_threshold)`.
#'   Defaults to `0.2`.
#' @param center Logical, whether to center the data before PCA (passed to `stats::prcomp`). Defaults to `TRUE`.
#' @param scale. Logical, whether to scale the data before PCA (passed to `stats::prcomp`). Defaults to `FALSE`.
#' @param run_idx Optional. Numeric vector indicating run membership for each
#'   timepoint. Currently not used but reserved for future extensions like
#'   run-wise PCA.
#'
#' @return A matrix containing the extracted temporal principal components
#'   (timepoints x `n_pcs`). Returns `NULL` if critical errors occur (e.g.,
#'   `n_pcs` is invalid or no valid voxels for PCA).
#'
#' @examples
#' \dontrun{
#' # Basic Example
#' n_time <- 200
#' n_voxels <- 500
#' n_task_regs <- 2
#' n_motion_regs <- 6
#'
#' Y_fmri_example <- matrix(rnorm(n_time * n_voxels), n_time, n_voxels)
#' X_base_example <- matrix(rnorm(n_time * (n_task_regs + n_motion_regs)), n_time,
#'                         n_task_regs + n_motion_regs)
#'
#' gdlite_pcs <- ndx_extract_gdlite_pcs(
#'   Y_fmri = Y_fmri_example,
#'   X_nuisance_base = X_base_example,
#'   n_pcs = 5
#' )
#'
#' if (!is.null(gdlite_pcs)) {
#'   print(paste("Extracted GD-Lite PCs dimensions:", dim(gdlite_pcs)[1], "x", dim(gdlite_pcs)[2]))
#' }
#'
#' # With a voxel mask
#' voxel_mask_example <- sample(c(TRUE, FALSE), n_voxels, replace = TRUE, prob = c(0.3, 0.7))
#' gdlite_pcs_masked <- ndx_extract_gdlite_pcs(
#'   Y_fmri = Y_fmri_example,
#'   X_nuisance_base = X_base_example,
#'   n_pcs = 3,
#'   voxel_mask = voxel_mask_example
#' )
#' if (!is.null(gdlite_pcs_masked)) {
#'   print(paste("Extracted masked GD-Lite PCs dimensions:",
#'                dim(gdlite_pcs_masked)[1], "x", dim(gdlite_pcs_masked)[2]))
#' }
#' }
#' @import stats
#' @import matrixStats
#' @export
ndx_extract_gdlite_pcs <- function(Y_fmri,
                                   X_nuisance_base,
                                   n_pcs,
                                   voxel_mask = NULL,
                                   min_variance_voxels = NULL,
                                   max_variance_voxels_pca_prop = 0.2,
                                   center = TRUE,
                                   scale. = FALSE,
                                   run_idx = NULL) {

  # --- 1. Input Validation ---
  if (!is.matrix(Y_fmri) || !is.numeric(Y_fmri)) {
    stop("`Y_fmri` must be a numeric matrix.")
  }
  if (nrow(Y_fmri) == 0 || ncol(Y_fmri) == 0) {
    stop("`Y_fmri` must not be empty.")
  }
  n_timepoints <- nrow(Y_fmri)
  n_total_voxels <- ncol(Y_fmri)

  if (!is.matrix(X_nuisance_base) || !is.numeric(X_nuisance_base)) {
    stop("`X_nuisance_base` must be a numeric matrix.")
  }
  if (nrow(X_nuisance_base) != n_timepoints) {
    stop("`X_nuisance_base` must have the same number of rows as `Y_fmri`.")
  }
  if (ncol(X_nuisance_base) == 0) {
    warning("`X_nuisance_base` has zero columns. Residuals will be `Y_fmri` itself (if demeaned).")
  }
  
  if (!is.numeric(n_pcs) || length(n_pcs) != 1 || n_pcs <= 0 || n_pcs >= n_timepoints || n_pcs >= n_total_voxels ) {
    warning(sprintf("Invalid `n_pcs` (%s). Must be > 0 and < min(timepoints, voxels_for_pca). Setting n_pcs to min(5, default_n_pcs).", as.character(n_pcs)))
    n_pcs <- min(5, floor(min(n_timepoints, n_total_voxels)/2), na.rm = TRUE)
    if (n_pcs <=0) {
        warning("Cannot determine a valid n_pcs. Returning NULL.")
        return(NULL)
    }
  }
  n_pcs <- as.integer(n_pcs)


  if (is.null(min_variance_voxels)) {
    min_variance_voxels <- max(100L, n_pcs * 5L)
  }
  if (!is.numeric(min_variance_voxels) || length(min_variance_voxels) != 1 || min_variance_voxels <= 0) {
    stop("`min_variance_voxels` must be a positive integer.")
  }
  min_variance_voxels <- as.integer(min_variance_voxels)

  if (!is.numeric(max_variance_voxels_pca_prop) || length(max_variance_voxels_pca_prop) != 1 ||
      max_variance_voxels_pca_prop <= 0 || max_variance_voxels_pca_prop > 1) {
    stop("`max_variance_voxels_pca_prop` must be a numeric value between 0 (exclusive) and 1 (inclusive).")
  }

  # --- 2. Calculate Residuals ---
  # Ensure X_nuisance_base has an intercept if not already present for lm.fit
  # However, lm.fit handles residuals correctly even without explicit intercept if X is full rank.
  # For safety and explicit control, one might consider adding an intercept if mean-centering of Y is desired
  # but not guaranteed by X_nuisance_base.
  # Here, we assume X_nuisance_base is appropriately specified by the user.
  
  if (ncol(X_nuisance_base) > 0) {
    residuals_for_pca <- tryCatch({
      stats::lm.fit(x = X_nuisance_base, y = Y_fmri)$residuals
    }, error = function(e) {
      warning(paste("Error in lm.fit for initial residual calculation:", e$message, ". Returning NULL."))
      return(NULL)
    })
    if (is.null(residuals_for_pca)) return(NULL)
  } else { # No nuisance regressors, PCA on (potentially centered/scaled) Y_fmri
    residuals_for_pca <- Y_fmri
  }
  
  if (any(!is.finite(residuals_for_pca))) {
      warning("Non-finite values found in residuals prior to PCA. Replacing with 0.")
      residuals_for_pca[!is.finite(residuals_for_pca)] <- 0
  }


  # --- 3. Select Voxels for PCA ---
  selected_voxel_indices <- NULL
  if (!is.null(voxel_mask)) {
    if (is.logical(voxel_mask)) {
      if (length(voxel_mask) != n_total_voxels) {
        stop("Logical `voxel_mask` length must match the number of voxels in `Y_fmri`.")
      }
      selected_voxel_indices <- which(voxel_mask)
    } else if (is.numeric(voxel_mask)) {
      if (any(voxel_mask < 1) || any(voxel_mask > n_total_voxels)) {
        stop("Numeric `voxel_mask` contains out-of-bounds indices.")
      }
      selected_voxel_indices <- unique(as.integer(voxel_mask))
    } else {
      stop("`voxel_mask` must be logical, numeric, or NULL.")
    }
    if (length(selected_voxel_indices) == 0) {
      warning("`voxel_mask` resulted in zero selected voxels. Returning NULL.")
      return(NULL)
    }
  } else {
    # Select a subset of high-variance voxels from residuals
    if (n_total_voxels <= min_variance_voxels) {
      selected_voxel_indices <- 1:n_total_voxels
    } else {
      voxel_variances <- matrixStats::colVars(residuals_for_pca, na.rm = TRUE)
      if (all(!is.finite(voxel_variances)) || all(voxel_variances <= 1e-9) ){
          warning("All voxel variances are zero or non-finite after residualization. Cannot select high-variance voxels. Using first min_variance_voxels if possible.")
          selected_voxel_indices <- 1:min(n_total_voxels, min_variance_voxels)
      } else {
          # Protect against all zero variances after some operations
          finite_variances <- is.finite(voxel_variances) & voxel_variances > 1e-9
          if(!any(finite_variances)) {
              warning("No voxels with finite, positive variance found. Using first min_variance_voxels if possible.")
              selected_voxel_indices <- 1:min(n_total_voxels, min_variance_voxels)
          } else {
              num_candidate_voxels <- max(min_variance_voxels, floor(n_total_voxels * max_variance_voxels_pca_prop))
              num_to_select <- min(n_total_voxels, num_candidate_voxels, sum(finite_variances))
              
              # Order by variance (descending) only among finite positive variance voxels
              ordered_indices <- order(voxel_variances[finite_variances], decreasing = TRUE)
              selected_voxel_indices <- which(finite_variances)[ordered_indices[1:num_to_select]]
          }
      }
    }
  }
  
  if (length(selected_voxel_indices) == 0) {
    warning("No voxels selected for PCA. Returning NULL.")
    return(NULL)
  }
  
  num_voxels_for_pca <- length(selected_voxel_indices)
  if (n_pcs >= num_voxels_for_pca) {
    warning(sprintf("`n_pcs` (%d) is >= number of selected voxels for PCA (%d). Adjusting n_pcs to %d.",
                    n_pcs, num_voxels_for_pca, max(1L, num_voxels_for_pca - 1L)))
    n_pcs <- max(1L, num_voxels_for_pca - 1L) # PCA rank is min(N,P)-1 if centered. If P < N, rank is P-1.
    if (n_pcs == 0 && num_voxels_for_pca > 0) n_pcs <- 1 # if only 1 voxel, can still extract 1 PC (itself)
    if (n_pcs == 0) {
        warning("Adjusted n_pcs is 0. Returning NULL.")
        return(NULL)
    }
  }
   if (n_pcs >= n_timepoints) {
    warning(sprintf("`n_pcs` (%d) is >= number of timepoints (%d). Adjusting n_pcs to %d.",
                    n_pcs, n_timepoints, max(1L, n_timepoints - 1L)))
    n_pcs <- max(1L, n_timepoints - 1L)
     if (n_pcs == 0) {
        warning("Adjusted n_pcs is 0 due to timepoints. Returning NULL.")
        return(NULL)
    }
  }


  data_for_pca <- residuals_for_pca[, selected_voxel_indices, drop = FALSE]

  # --- 4. Perform PCA ---
  pca_results <- tryCatch({
    stats::prcomp(data_for_pca, rank. = n_pcs, center = center, scale. = scale.)
  }, error = function(e) {
    warning(paste("Error during PCA execution (stats::prcomp):", e$message, ". Returning NULL."))
    return(NULL)
  })

  if (is.null(pca_results) || is.null(pca_results$x)) {
    warning("PCA did not return valid components. Returning NULL.")
    return(NULL)
  }

  # --- 5. Extract and Return PCs ---
  # prcomp returns PCs in $x
  extracted_pcs <- pca_results$x[, 1:min(n_pcs, ncol(pca_results$x)), drop = FALSE]
  
  if (ncol(extracted_pcs) != n_pcs) {
      warning(sprintf("PCA returned %d components, but %d were requested. This might be due to data rank.", ncol(extracted_pcs), n_pcs))
  }

  return(extracted_pcs)
}

#' @title Calculate Leave-One-Run-Out Cross-Validated R-squared
#'
#' @description
#' Computes the R-squared for each voxel using leave-one-run-out (LORO)
#' cross-validation. For each run, a GLM is fit to all other (training) runs,
#' and then used to predict the data in the left-out (test) run. The R-squared
#' is then calculated for each voxel based on the residuals from this prediction
#' in the test run.
#'
#' @param Y_fmri Numeric matrix of fMRI data (timepoints x voxels), concatenated
#'   across all runs.
#' @param X_design Numeric matrix, the design matrix (timepoints x regressors) for the GLM.
#'   This should correspond to the full `Y_fmri` matrix (concatenated across runs).
#' @param run_idx Numeric vector indicating run membership for each timepoint (row)
#'   in `Y_fmri` and `X_design`.
#'
#' @return A numeric vector of length `ncol(Y_fmri)` containing the LORO
#'   cross-validated R-squared value for each voxel. If a voxel has zero variance
#'   in a test run, or if the GLM fit fails for a fold, its R-squared might be 0 or NA.
#'
#' @importFrom stats lm.fit var
#' @keywords internal
#' @export
cv_r2_loro <- function(Y_fmri, X_design, run_idx) {
  if (!is.matrix(Y_fmri) || !is.numeric(Y_fmri)) stop("`Y_fmri` must be a numeric matrix.")
  if (!is.matrix(X_design) || !is.numeric(X_design)) stop("`X_design` must be a numeric matrix.")
  if (!is.numeric(run_idx)) stop("`run_idx` must be a numeric vector.")

  if (nrow(Y_fmri) != nrow(X_design)) {
    stop("Number of rows in `Y_fmri` and `X_design` must match.")
  }
  if (nrow(Y_fmri) != length(run_idx)) {
    stop("Length of `run_idx` must match the number of rows in `Y_fmri`.")
  }
  if (ncol(X_design) == 0) {
    warning("`X_design` has zero columns. LORO R2 cannot be meaningfully computed. Returning zeros.")
    return(rep(0, ncol(Y_fmri)))
  }

  unique_runs <- sort(unique(run_idx))
  n_runs <- length(unique_runs)
  n_voxels <- ncol(Y_fmri)
  n_timepoints_total <- nrow(Y_fmri)

  if (n_runs <= 1) {
    warning("LORO R2 requires at least 2 runs. Cannot perform cross-validation. Returning NAs.")
    return(rep(NA_real_, n_voxels))
  }

  # Initialize matrix to store squared errors for each voxel from its test run prediction
  sum_sq_residuals_test <- numeric(n_voxels)
  sum_sq_total_test <- numeric(n_voxels)
  # Keep track of voxels that have valid calculations (to handle all-NA cases)
  voxel_is_valid <- logical(n_voxels)

  for (r_idx in seq_along(unique_runs)) {
    current_test_run <- unique_runs[r_idx]
    if (getOption("fmridenoise.verbose_gdlite", FALSE)) {
        message(sprintf("  cv_r2_loro: Processing fold %d/%d (leaving out run %s)", r_idx, n_runs, as.character(current_test_run)))
    }

    train_mask <- run_idx != current_test_run
    test_mask <- run_idx == current_test_run

    if (sum(train_mask) == 0 || sum(test_mask) == 0) {
      warning(sprintf("Skipping run %s in LORO CV: no train or test samples.", as.character(current_test_run)))
      next
    }

    Y_train <- Y_fmri[train_mask, , drop = FALSE]
    X_train <- X_design[train_mask, , drop = FALSE]
    Y_test <- Y_fmri[test_mask, , drop = FALSE]
    X_test <- X_design[test_mask, , drop = FALSE]

    # Fit GLM on training data
    # Check for rank deficiency in X_train for this fold
    qr_X_train <- qr(X_train)
    if (qr_X_train$rank < ncol(X_train)) {
        if (getOption("fmridenoise.verbose_gdlite", FALSE)) {
            message(sprintf("    Rank deficient X_train for fold %d (run %s out). Rank = %d, Cols = %d", 
                            r_idx, as.character(current_test_run), qr_X_train$rank, ncol(X_train)))
        }
        # lm.fit should still proceed with a pseudo-inverse
    }
    
    betas_train <- tryCatch({
      stats::lm.fit(x = X_train, y = Y_train)$coefficients
    }, error = function(e) {
      warning(sprintf("lm.fit failed for training data in LORO fold (run %s out): %s", as.character(current_test_run), e$message))
      return(NULL)
    })

    if (is.null(betas_train)) {
      # Cannot proceed with this fold for prediction if beta estimation failed.
      # Voxels in this test run won't get updated R2 from this fold.
      next
    }
    
    # Handle cases where some voxels might have all NA betas if Y_train columns were all NA
    # This can happen if a voxel is entirely outside a mask for all training runs.
    # Such voxels won't contribute to R2 calculation for this fold.
    # `betas_train` will have NA columns for such voxels if `lm.fit` encounters all NA y.

    # Predict on test data
    Y_pred_test <- X_test %*% betas_train
    residuals_test <- Y_test - Y_pred_test

    # Calculate Sum of Squared Residuals (RSS) and Total Sum of Squares (TSS) for each voxel in the test set
    for (v_idx in 1:n_voxels) {
      y_v_test <- Y_test[, v_idx]
      res_v_test <- residuals_test[, v_idx]
      
      # Ensure there are valid (non-NA) observations for this voxel in this test run
      valid_obs_mask_v_test <- !is.na(y_v_test)
      if (sum(valid_obs_mask_v_test) < 2) { # Need at least 2 points to calculate variance for TSS
          next # Skip this voxel for this fold if not enough data
      }
      
      y_v_test_valid <- y_v_test[valid_obs_mask_v_test]
      res_v_test_valid <- res_v_test[valid_obs_mask_v_test]
      
      # If after NA removal for y, residuals are all NA (e.g. because betas_train[,v_idx] was NA)
      if(all(is.na(res_v_test_valid))) {
          next
      }

      current_rss_v_test <- sum(res_v_test_valid^2, na.rm = TRUE)
      mean_y_v_test_valid <- mean(y_v_test_valid, na.rm = TRUE)
      current_tss_v_test <- sum((y_v_test_valid - mean_y_v_test_valid)^2, na.rm = TRUE)
      
      # Accumulate (this is incorrect for LORO, R2 is per-fold, then assigned)
      # For LORO, R2 for a voxel is computed *when its run is the test set*.
      # So, we don't accumulate. We store these per-voxel for this specific left-out run.
      # The final r2_cv_voxels will be populated such that r2_cv_voxels[v] is the R2 calculated when
      # the run(s) containing voxel v data was part of Y_test.
      # This needs rethinking: cv_r2_loro should return one R2 per voxel.
      # The R2 for voxel `v` is the R2 obtained when *its own run* was the test data.
      # The loop structure implies we are calculating R2 for ALL voxels using the current fold.
      # This is fine, but how to store it? The R2 for voxel `v` belonging to run `R` is the one
      # calculated when run `R` is `current_test_run`.
      
      # Let's re-initialize r2_cv_voxels inside the loop and assign to relevant voxels.
      # This is still not quite right. The output must be one R2 per voxel.
      # The standard way: for each voxel, its LORO R2 is the R2 obtained from the fold where its run was held out.
      # Thus, after this loop, `r2_cv_values_for_voxels_in_test_run` is what we need for this subset of voxels.
      
      # Store for later: total sum of squares and residual sum of squares for *this specific test run data*
      # This is what the GLMdenoise paper implies: R^2(v) = 1 - RSS_test(v) / TSS_test(v)
      # where test is when run(v) is left out.
      
      # If current_tss_v_test is near zero, R2 is 0 (or undefined, treat as 0 for robustness)
      if (abs(current_tss_v_test) < 1e-9) {
        # R2 for this voxel in this fold is 0
      } else {
        # R2 for this voxel, for this fold (when its run was left out)
        # sum_sq_residuals_test[v_idx] and sum_sq_total_test[v_idx] will be populated *only* for voxels in current_test_run
        # This means we need to know which voxels belong to which runs if runs have different voxels (not typical)
        # Assuming all voxels are present in all runs for now.
        # The current loop calculates R2 for ALL voxels based on the test data of the CURRENT FOLD.
        # This means for voxels in run `q` their `r2_cv` is taken from the fold where run `q` was left out.
        
        # This implies that for all voxels `v_idx`, when their run is the `current_test_run`,
        # their `sum_sq_residuals_test[v_idx]` and `sum_sq_total_test[v_idx]` should be updated.
        
        # If a voxel `v_idx` is part of `current_test_run` (which all are, conceptually, as Y_test is [time_in_test_run x all_voxels])
        # then update its specific RSS and TSS based on this fold where its run was tested.
        sum_sq_residuals_test[v_idx] <- sum_sq_residuals_test[v_idx] + current_rss_v_test
        sum_sq_total_test[v_idx] <- sum_sq_total_test[v_idx] + current_tss_v_test
        voxel_is_valid[v_idx] <- TRUE # Mark that this voxel had at least one valid LORO calculation
      }
    }
  } # End LORO run loop

  # Calculate final R2_cv for each voxel
  r2_cv_voxels <- numeric(n_voxels)
  for (v_idx in 1:n_voxels) {
    if (!voxel_is_valid[v_idx] || abs(sum_sq_total_test[v_idx]) < 1e-9 * n_runs) { # Check if TSS accumulated is too small
      r2_cv_voxels[v_idx] <- 0 # Or NA_real_ if no valid LORO calculations or TSS is zero across all folds relevant to it
    } else {
      # The sums are over predictions of ALL voxels in EACH test run.
      # For GLMdenoise, R^2_CV(v) = 1 - sum(Residuals_test^2(v)) / sum((Data_test(v) - mean(Data_test(v)))^2)
      # where the sum is over all TRs *when voxel v was in a test set*. This is what the accumulation does.
      r2_val <- 1 - (sum_sq_residuals_test[v_idx] / sum_sq_total_test[v_idx])
      # Clamp R2 values: can be < 0 if model predicts worse than mean. GLMdenoise seems to cap at 0.
      if (!is.finite(r2_val) || r2_val < 0) r2_val <- 0
      r2_cv_voxels[v_idx] <- r2_val
    }
  }
  
  # Ensure R2 is not slightly negative due to floating point issues if model is very poor
  r2_cv_voxels[r2_cv_voxels < 0] <- 0
  return(r2_cv_voxels)
}

#' @title Calculate Temporal Signal-to-Noise Ratio (tSNR)
#'
#' @description
#' Computes the temporal Signal-to-Noise Ratio (tSNR) for each voxel in an fMRI dataset.
#' tSNR is typically calculated as the mean of a voxel's time series divided by its
#' standard deviation over time.
#'
#' @param Y_data Numeric matrix of fMRI data or residuals (timepoints x voxels).
#' @param detrend Logical, if TRUE, linearly detrend each voxel's time series before
#'   calculating mean and standard deviation. Defaults to `FALSE` as input data
#'   (e.g. residuals from a GLM that included trends) might already be detrended.
#' @param robust Logical, if `TRUE`, use median for mean and MAD for standard deviation.
#'   Default is `FALSE`.
#'
#' @return A numeric vector of length `ncol(Y_data)` containing the tSNR value for
#'   each voxel. Returns `NA` for voxels with zero standard deviation or if issues occur.
#'
#' @importFrom matrixStats colMeans2 colSds colMedians colMads
#' @keywords internal
#' @export
calculate_tsnr <- function(Y_data, detrend = FALSE, robust = FALSE) {
  if (!is.matrix(Y_data) || !is.numeric(Y_data)) {
    stop("`Y_data` must be a numeric matrix.")
  }
  if (ncol(Y_data) == 0 || nrow(Y_data) < 2) {
    warning("`Y_data` has zero voxels or insufficient timepoints (<2). Returning empty numeric.")
    return(numeric(0))
  }

  n_timepoints <- nrow(Y_data)
  time_vector <- seq_len(n_timepoints)
  residuals <- Y_data

  if (detrend) {
    t_mean <- mean(time_vector)
    time_centered <- time_vector - t_mean

    mean_y <- matrixStats::colMeans2(Y_data, na.rm = TRUE)
    y_centered <- sweep(Y_data, 2, mean_y, "-")
    valid_mask <- !is.na(Y_data)

    numer <- colSums(time_centered * y_centered, na.rm = TRUE)
    denom <- colSums((time_centered^2) * valid_mask, na.rm = TRUE)
    slope <- numer / denom
    slope[!is.finite(slope)] <- NA_real_
    intercept <- mean_y - slope * t_mean

    pred <- outer(time_vector, slope)
    pred <- sweep(pred, 2, intercept, "+")
    residuals <- Y_data - pred
  }

  valid_counts <- colSums(!is.na(residuals))

  if (robust) {
    mean_vals <- matrixStats::colMedians(residuals, na.rm = TRUE)
    sd_vals <- matrixStats::colMads(residuals, na.rm = TRUE)
  } else {
    mean_vals <- matrixStats::colMeans2(residuals, na.rm = TRUE)
    sd_vals <- matrixStats::colSds(residuals, na.rm = TRUE)
  }

  sd_vals[valid_counts < 2 | sd_vals < 1e-9 | !is.finite(sd_vals)] <- NA_real_

  tsnr_values <- mean_vals / sd_vals
  tsnr_values[!is.finite(tsnr_values)] <- NA_real_

  return(tsnr_values)
}

#' @title Select Optimal Number of Principal Components for GLMdenoise
#'
#' @description
#' This function determines the optimal number of GLMdenoise principal components (PCs)
#' to include in a GLM. It iterates from 1 to `K_max` PCs, appends them to a
#' base design matrix (`X_base_design`), and recalculates leave-one-run-out (LORO)
#' cross-validated R-squared (`r2_cv`). The optimal number of PCs (`K_star`) is
#' chosen as the one that maximizes the median `r2_cv` across a predefined set of
#' "good voxels".
#'
#' @param Y_fmri Numeric matrix of fMRI data (timepoints x voxels), concatenated across runs.
#' @param X_base_design Numeric matrix, the base design matrix (e.g., task regressors,
#'   motion parameters, polynomial trends) to which PCs will be appended.
#'   Must have the same number of rows as `Y_fmri`.
#' @param all_pcs Numeric matrix (timepoints x `K_max`), containing all candidate PCs
#'   extracted from noise-pool residuals (e.g., output of `ndx_extract_gdlite_pcs`
#'   called with `n_pcs = K_max`).
#' @param run_idx Numeric vector indicating run membership for each timepoint in `Y_fmri`.
#' @param good_voxel_mask Logical vector of length `ncol(Y_fmri)`, where TRUE indicates
#'   a "good voxel" over which the median LORO R-squared will be calculated to select `K_star`.
#'   These are typically voxels with a reasonable initial R-squared from a simpler model.
#' @param K_max Integer, the maximum number of PCs to consider. Should match `ncol(all_pcs)`.
#'   If `NULL`, it's derived from `ncol(all_pcs)`.
#' @param min_K Integer, the minimum number of PCs to consider in the loop. Default is 0 (model with no PCs).
#'
#' @return A list containing:
#'   - `K_star`: The optimal number of PCs selected.
#'   - `selected_pcs`: The matrix of selected PCs (`all_pcs[, 1:K_star, drop = FALSE]`).
#'     If `K_star` is 0, this will be `NULL` or an empty matrix.
#'   - `median_r2cv_values`: A numeric vector of median R-squared_cv values for each `K` tested.
#'   - `all_r2cv_per_K`: A list where each element is the full r2_cv vector for all voxels for a given K.
#'
#' @importFrom stats median
#' @keywords internal
#' @export
select_optimal_k_gdlite <- function(Y_fmri,
                                      X_base_design,
                                      all_pcs,
                                      run_idx,
                                      good_voxel_mask,
                                      K_max = NULL,
                                      min_K = 0) {

  # --- Input Validation ---
  if (!is.matrix(Y_fmri) || !is.numeric(Y_fmri)) stop("`Y_fmri` must be a numeric matrix.")
  if (!is.matrix(X_base_design) || !is.numeric(X_base_design)) stop("`X_base_design` must be a numeric matrix.")
  if (!is.matrix(all_pcs) || !is.numeric(all_pcs)) stop("`all_pcs` must be a numeric matrix.")
  if (!is.numeric(run_idx)) stop("`run_idx` must be a numeric vector.")
  if (!is.logical(good_voxel_mask) || length(good_voxel_mask) != ncol(Y_fmri)) {
    stop("`good_voxel_mask` must be a logical vector matching the number of voxels in `Y_fmri`.")
  }
  if (sum(good_voxel_mask) == 0) {
    warning("`good_voxel_mask` has no TRUE values. Cannot select optimal K. Returning K_star = 0.")
    return(list(K_star = 0, selected_pcs = NULL, median_r2cv_values = numeric(0), all_r2cv_per_K = list()))
  }

  if (nrow(Y_fmri) != nrow(X_base_design) || nrow(Y_fmri) != nrow(all_pcs)) {
    stop("`Y_fmri`, `X_base_design`, and `all_pcs` must have the same number of rows (timepoints).")
  }
  if (nrow(Y_fmri) != length(run_idx)) {
    stop("Length of `run_idx` must match the number of rows in `Y_fmri`.")
  }
  if (is.null(K_max)) {
    K_max <- ncol(all_pcs)
  }
  if (!is.numeric(K_max) || K_max < 0 || K_max > ncol(all_pcs)) {
    stop(sprintf("Invalid `K_max` (%s). Must be between 0 and ncol(all_pcs) (%d).", 
                 as.character(K_max), ncol(all_pcs)))
  }
  if (!is.numeric(min_K) || min_K < 0 || min_K > K_max) {
      stop(sprintf("Invalid `min_K` (%s). Must be between 0 and K_max (%d).", as.character(min_K), K_max))
  }
  K_max <- as.integer(K_max)
  min_K <- as.integer(min_K)
  
  num_loops <- K_max - min_K + 1
  median_r2cv_values <- numeric(num_loops)
  all_r2cv_per_K_list <- vector("list", length = num_loops)
  k_values_tested <- min_K:K_max
  loop_idx <- 0

  if (getOption("fmridenoise.verbose_gdlite", FALSE)) {
      message(sprintf("  select_optimal_k_gdlite: Testing K from %d to %d.", min_K, K_max))
  }

  for (k_current_pcs in k_values_tested) {
    loop_idx <- loop_idx + 1
    if (getOption("fmridenoise.verbose_gdlite", FALSE)) {
      message(sprintf("    Testing K = %d principal components...", k_current_pcs))
    }

    X_current_design <- X_base_design
    if (k_current_pcs > 0) {
      current_pcs_to_add <- all_pcs[, 1:k_current_pcs, drop = FALSE]
      # Ensure PCs have unique names to avoid issues if X_base_design has similar names
      colnames(current_pcs_to_add) <- paste0("GDLitePC_selK_", seq_len(ncol(current_pcs_to_add)))
      X_current_design <- cbind(X_base_design, current_pcs_to_add)
      # It might be good to ensure X_current_design colnames are globally unique if X_base_design is complex
      if (any(duplicated(colnames(X_current_design)))) {
          colnames(X_current_design) <- make.names(colnames(X_current_design), unique=TRUE)
      }
    }
    
    if (ncol(X_current_design) >= nrow(Y_fmri)) {
        warning(sprintf("K = %d: Number of regressors (%d) >= number of timepoints (%d). This may lead to unstable LORO R2. Skipping.", 
                        k_current_pcs, ncol(X_current_design), nrow(Y_fmri)))
        median_r2cv_values[loop_idx] <- -Inf # Or NA, to ensure it's not picked
        all_r2cv_per_K_list[[loop_idx]] <- rep(NA_real_, ncol(Y_fmri))
        next
    }

    r2_cv_current_k <- cv_r2_loro(Y_fmri, X_current_design, run_idx)
    all_r2cv_per_K_list[[loop_idx]] <- r2_cv_current_k
    
    if (any(!is.na(r2_cv_current_k[good_voxel_mask]))){
        median_r2cv_values[loop_idx] <- stats::median(r2_cv_current_k[good_voxel_mask], na.rm = TRUE)
    } else {
        # If all R2cv values in good voxels are NA for this K, assign a very low value
        median_r2cv_values[loop_idx] <- -Inf 
    }
    
    if (getOption("fmridenoise.verbose_gdlite", FALSE)) {
      message(sprintf("      K = %d, Median R2cv (good voxels): %.4f", k_current_pcs, median_r2cv_values[loop_idx]))
    }
  }

  if (all(is.infinite(median_r2cv_values)) || all(is.na(median_r2cv_values))) {
      warning("All median R2cv values were -Inf or NA. Cannot determine optimal K using R2cv.")
      # Fallback for K_star selection when R2cv is not informative (e.g. single run)
      if (getOption("fmridenoise.verbose_gdlite", FALSE)) {
          message(sprintf("  Attempting fallback K_star selection. k_values_tested range: %d to %d. Actual extracted PCs: %d", 
                          min_K, K_max, if(is.null(all_pcs)) 0 else ncol(all_pcs)))
      }
      if (!is.null(all_pcs) && ncol(all_pcs) > 0) {
          # If min_K_optimal_selection is > 0 and <= available PCs, use it.
          # Otherwise, pick a small number, e.g., 1, if available.
          if (min_K > 0 && min_K <= ncol(all_pcs)) {
              K_star <- min_K
              if (getOption("fmridenoise.verbose_gdlite", FALSE)) message(sprintf("  Fallback K_star: Set to min_K = %d", K_star))
          } else if (ncol(all_pcs) >= 1) {
              K_star <- 1L # Fallback to 1 PC if available and min_K was 0 or too high
              if (getOption("fmridenoise.verbose_gdlite", FALSE)) message(sprintf("  Fallback K_star: Set to 1 (min_K was %d)", min_K))
          } else {
              K_star <- 0 # No PCs available
              if (getOption("fmridenoise.verbose_gdlite", FALSE)) message("  Fallback K_star: No candidate PCs available, set to 0.")
          }
          # Ensure K_star does not exceed K_max (actual number of PCs extracted)
          if (!is.null(all_pcs)) K_star <- min(K_star, ncol(all_pcs))
          
      } else {
          K_star <- 0 # No candidate PCs to select from
          if (getOption("fmridenoise.verbose_gdlite", FALSE)) message("  Fallback K_star: No candidate PCs provided, set to 0.")
      }
      if (K_star == 0) warning("Fallback K_star selection resulted in K_star = 0.")
      
  } else {
      best_loop_idx <- which.max(median_r2cv_values)
      if (length(best_loop_idx) > 1) best_loop_idx <- best_loop_idx[1] # Pick first if multiple maxima
      K_star <- k_values_tested[best_loop_idx]
  }
  
  if (getOption("fmridenoise.verbose_gdlite", FALSE)) {
      message(sprintf("  Optimal K_star selected: %d (Median R2cv: %.4f)", K_star, median_r2cv_values[which(k_values_tested == K_star)]))
  }

  selected_pcs_out <- NULL
  if (K_star > 0) {
    selected_pcs_out <- all_pcs[, 1:K_star, drop = FALSE]
    colnames(selected_pcs_out) <- paste0("GDLitePC_final_", seq_len(ncol(selected_pcs_out))) # Consistent final naming
  }
  
  names(all_r2cv_per_K_list) <- paste0("K=", k_values_tested)
  names(median_r2cv_values) <- paste0("K=", k_values_tested)

  return(list(
    K_star = K_star,
    selected_pcs = selected_pcs_out,
    median_r2cv_values = median_r2cv_values,
    all_r2cv_per_K = all_r2cv_per_K_list
  ))
}

#' @title Run a GLMdenoise-Lite Analysis
#'
#' @description
#' Implements a simplified version of the GLMdenoise procedure to identify
#' data-driven noise regressors (principal components from residuals).
#' The steps include:
#' 1. Constructing a base GLM design matrix (`X_gd`) including task regressors,
#'    polynomial trends, and motion parameters.
#' 2. Calculating leave-one-run-out cross-validated R-squared (`r2_cv_initial`)
#'    based on `X_gd`.
#' 3. Calculating temporal Signal-to-Noise Ratio (tSNR) for each voxel.
#' 4. Defining a "noise pool" of voxels based on low `r2_cv_initial` and adequate tSNR.
#' 5. Extracting candidate principal components (`pcs_all`) from the residuals of
#'    `Y_fmri ~ X_gd` within the noise pool.
#' 6. Defining a "good voxel pool" (e.g., `r2_cv_initial` > threshold) for PC selection.
#' 7. Selecting the optimal number of PCs (`K_star`) by maximizing the median `r2_cv`
#'    (from models `Y_fmri ~ X_gd + K_selected_PCs`) within the good voxel pool.
#' 8. Optionally, performing a final GLM fit using `X_gd` and the `selected_pcs`.
#'
#' @param Y_fmri Numeric matrix of fMRI data (timepoints x voxels), concatenated across runs.
#' @param events A data frame describing experimental events, required to build task regressors.
#'   Must contain columns `onsets`, `durations`, `condition`, `blockids`.
#' @param run_idx Numeric vector indicating run membership for each timepoint.
#' @param TR Numeric, repetition time in seconds.
#' @param motion_params Optional numeric matrix of motion parameters (timepoints x N).
#'   If provided, these will be included in `X_gd`.
#' @param poly_degree Integer or NULL. Degree of polynomial trends to include per run.
#'   If NULL, no polynomial trends are added. From `fmrireg` conventions, 0 means run-specific intercepts.
#'   Default is 1 (linear trend per run).
#' @param k_max Integer, maximum number of PCs to extract and consider. Default is 30.
#' @param r2_thresh_noise_pool Numeric, R-squared_cv threshold to define noise pool voxels
#'   (voxels with `r2_cv_initial < r2_thresh_noise_pool` are part of noise pool if tSNR is also met).
#'   Default is 0.05.
#' @param tsnr_thresh_noise_pool Numeric, tSNR threshold to define noise pool voxels
#'   (voxels with `tsnr > tsnr_thresh_noise_pool` are part of noise pool if R2_cv is also met).
#'   Default is 30.
#' @param r2_thresh_good_voxels Numeric, R-squared_cv threshold to define the "good voxel pool"
#'   used for selecting `K_star` (median R2_cv is maximized over voxels with
#'   `r2_cv_initial >= r2_thresh_good_voxels`). Default is 0.05.
#' @param min_K_optimal_selection Integer, minimum number of PCs to consider when selecting optimal K.
#'   Default is 0 (allows model with no PCs to be chosen).
#' @param detrend_for_tsnr Logical, whether to detrend data before tSNR calculation. Default is TRUE.
#' @param perform_final_glm Logical, if TRUE, a final GLM is fit using `X_gd` and the
#'   selected PCs, and its coefficients and residuals are returned. Default is TRUE.
#'
#' @return A list containing various outputs from the GLMdenoise-Lite procedure:
#'   - `selected_pcs`: The matrix of `K_star` selected principal components.
#'   - `K_star`: The optimal number of PCs selected.
#'   - `X_gd`: The base design matrix used (task, motion, polynomials).
#'   - `r2_cv_initial`: The initial leave-one-run-out R-squared values for `X_gd`.
#'   - `tsnr_values`: The tSNR value for each voxel.
#'   - `noise_pool_mask`: Logical vector indicating voxels selected for the noise pool.
#'   - `good_voxel_mask_for_k_selection`: Logical vector indicating voxels used for K selection.
#'   - `median_r2cv_by_K`: Median R-squared_cv values for each number of PCs tested.
#'   - `all_r2cv_per_K`: List of all R-squared_cv vectors for each K.
#'   - `betas_gdlite`: Coefficients from the final GLM (if `perform_final_glm = TRUE`).
#'   - `residuals_gdlite`: Residuals from the final GLM (if `perform_final_glm = TRUE`).
#'   - `all_candidate_pcs`: The full set of `k_max` PCs extracted from the noise pool.
#'
#' @importFrom fmrireg hrf_spmg1 event_model design_matrix sampling_frame baseline_model
#' @importFrom stats lm.fit
#' @export
ndx_run_gdlite <- function(Y_fmri,
                             events,
                             run_idx,
                             TR,
                             motion_params = NULL,
                             poly_degree = 1, # Consistent with fmrireg::baseline_model, 0=intercept, 1=linear, etc.
                             k_max = 30L,
                             r2_thresh_noise_pool = 0.05,
                             tsnr_thresh_noise_pool = 30,
                             r2_thresh_good_voxels = 0.05,
                             min_K_optimal_selection = 0L,
                             detrend_for_tsnr = TRUE,
                             perform_final_glm = TRUE) {

  if (getOption("fmridenoise.verbose_gdlite", FALSE)) {
    message("Starting GLMdenoise-Lite Analysis...")
  }
  n_timepoints <- nrow(Y_fmri)
  n_voxels <- ncol(Y_fmri)

  blocks_events <- sort(unique(events$blockids))
  blocks_runs <- sort(unique(run_idx))
  if (!identical(blocks_events, blocks_runs)) {
    stop(sprintf(
      "events$blockids %s do not match run_idx %s",
      paste(blocks_events, collapse = ","),
      paste(blocks_runs, collapse = ",")
    ))
  }
  events <- events[order(factor(events$blockids, levels = unique(run_idx))), , drop = FALSE]

  # --- Step 1: Build X_gd (Base Design Matrix) ---
  if (getOption("fmridenoise.verbose_gdlite", FALSE)) message("Step 1: Building base design matrix X_gd...")
  sf <- fmrireg::sampling_frame(blocklens = as.numeric(table(run_idx)), TR = TR)
  
  # Task model using canonical HRF (Glover + derivatives often used, spmg1 is simpler starter)
  # GLMdenoise typically uses Glover + time derivative + dispersion derivative.
  # For simplicity, starting with spmg1 (like ndx_initial_glm).
  # This could be a parameter: hrf_basis = "spmg1" or "glover_derivs"
  task_model_formula <- stats::as.formula("onsets ~ fmrireg::hrf(condition, basis='spmg1')")

  event_des <- fmrireg::event_model(task_model_formula,
                                    data = events,
                                    block = events$blockids,
                                    sampling_frame = sf)
  X_task <- fmrireg::design_matrix(event_des)
  if (is.null(X_task)) X_task <- matrix(0, nrow=n_timepoints, ncol=0)
  
  # Baseline model (polynomials per run, and run-specific intercepts if poly_degree allows)
  # fmrireg::baseline_model with basis="poly" and degree=poly_degree handles run intercepts correctly.
  # No nuisance regressors are included at this stage, so rely on the
  # default `nuisance_list = NULL` behavior. Passing an empty list would
  # cause an error inside `fmrireg::baseline_model` because it expects the
  # list length and row counts to match the sampling frame.
  baseline_des <- fmrireg::baseline_model(sframe = sf,
                                          basis = "poly",
                                          degree = poly_degree)
  X_poly_intercepts <- fmrireg::design_matrix(baseline_des)
  if (is.null(X_poly_intercepts)) X_poly_intercepts <- matrix(1, nrow=n_timepoints, ncol=1)

  X_gd_parts <- list(task = X_task, poly_intercepts = X_poly_intercepts)
  if (!is.null(motion_params)) {
    if (nrow(motion_params) == n_timepoints) {
      X_gd_parts$motion <- motion_params
    } else {
      warning("Motion parameters provided but row count does not match Y_fmri. Ignoring motion_params.")
    }
  }
  X_gd <- do.call(cbind, X_gd_parts)
  X_gd <- as.matrix(X_gd)
  storage.mode(X_gd) <- "numeric"
  if (any(duplicated(colnames(X_gd)))) {
      colnames(X_gd) <- make.names(colnames(X_gd), unique=TRUE)
  }
  if (ncol(X_gd) == 0) {
      warning("X_gd is empty. Cannot proceed with GLMdenoise-Lite.")
      return(NULL) # Or a more informative error list
  }
  if (getOption("fmridenoise.verbose_gdlite", FALSE)) message(sprintf("  X_gd built with %d regressors.", ncol(X_gd)))

  # --- Step 2: Initial LORO R-squared --- 
  if (getOption("fmridenoise.verbose_gdlite", FALSE)) message("Step 2: Calculating initial LORO R-squared (r2_cv_initial)...")
  r2_cv_initial <- cv_r2_loro(Y_fmri, X_gd, run_idx)
  is_single_run_or_cv_failed <- all(is.na(r2_cv_initial)) # Check if LORO CV returned all NAs

  # --- Step 3: Calculate tSNR --- 
  if (getOption("fmridenoise.verbose_gdlite", FALSE)) message("Step 3: Calculating tSNR...")
  tsnr_values <- calculate_tsnr(Y_fmri, detrend = detrend_for_tsnr)
  tsnr_values[is.na(tsnr_values)] <- 0 # Replace NA tSNR with 0 for mask creation

  # --- Step 4: Define Noise Pool Mask --- 
  if (getOption("fmridenoise.verbose_gdlite", FALSE)) message("Step 4: Defining noise pool mask...")
  if (is_single_run_or_cv_failed) {
    if (getOption("fmridenoise.verbose_gdlite", FALSE)) message("  Single run or LORO CV failed, using only TSNR for noise pool definition.")
    noise_pool_mask <- (tsnr_values > tsnr_thresh_noise_pool)
  } else {
    noise_pool_mask <- (r2_cv_initial < r2_thresh_noise_pool) & (tsnr_values > tsnr_thresh_noise_pool)
  }
  noise_pool_mask[is.na(noise_pool_mask)] <- FALSE # Ensure no NAs in mask
  if (sum(noise_pool_mask) == 0) {
    warning("Noise pool mask is empty (0 voxels). PCA will not be performed. Returning K_star=0.")
    # Fallback: return list indicating no PCs selected
    return(list(selected_pcs = NULL, K_star = 0, X_gd = X_gd, r2_cv_initial = r2_cv_initial,
                tsnr_values = tsnr_values, noise_pool_mask = noise_pool_mask, 
                good_voxel_mask_for_k_selection = rep(FALSE, n_voxels),
                median_r2cv_by_K = numeric(0), all_r2cv_per_K = list(),
                betas_gdlite = if(perform_final_glm) stats::lm.fit(X_gd, Y_fmri)$coefficients else NULL,
                residuals_gdlite = if(perform_final_glm) stats::lm.fit(X_gd, Y_fmri)$residuals else NULL,
                all_candidate_pcs = NULL))
  }
  if (getOption("fmridenoise.verbose_gdlite", FALSE)) message(sprintf("  Noise pool defined with %d voxels.", sum(noise_pool_mask)))

  # --- Step 5: Extract Candidate PCs from Noise Pool Residuals --- 
  if (getOption("fmridenoise.verbose_gdlite", FALSE)) message("Step 5: Extracting candidate PCs from noise pool residuals...")
  residuals_from_X_gd <- tryCatch({
    stats::lm.fit(x = X_gd, y = Y_fmri)$residuals
  }, error = function(e) {
    warning(paste("Error fitting GLM (X_gd on Y_fmri) to get residuals for PCA:", e$message))
    return(NULL)
  })
  if (is.null(residuals_from_X_gd)) {
      warning("Could not obtain residuals for PCA. Cannot proceed with PC extraction.")
      return(list(selected_pcs = NULL, K_star = 0, X_gd = X_gd, r2_cv_initial = r2_cv_initial,
                tsnr_values = tsnr_values, noise_pool_mask = noise_pool_mask, 
                good_voxel_mask_for_k_selection = rep(FALSE, n_voxels), # Placeholder
                median_r2cv_by_K = numeric(0), all_r2cv_per_K = list(),
                betas_gdlite = if(perform_final_glm) stats::lm.fit(X_gd, Y_fmri)$coefficients else NULL,
                residuals_gdlite = if(perform_final_glm) stats::lm.fit(X_gd, Y_fmri)$residuals else NULL,
                all_candidate_pcs = NULL))
  }
  
  # PCA is on residuals from X_gd, so X_nuisance_base for ndx_extract_gdlite_pcs should be empty or intercept-only.
  # Using an empty matrix for X_nuisance_base as residuals_from_X_gd are already residualized.
  all_candidate_pcs <- ndx_extract_gdlite_pcs(
    Y_fmri = residuals_from_X_gd, # PCA on residuals
    X_nuisance_base = matrix(0, nrow = n_timepoints, ncol = 0), # No further residualization needed
    n_pcs = k_max,
    voxel_mask = noise_pool_mask,
    center = TRUE, # Standard for PCA on residuals
    scale. = FALSE 
  )
  if (is.null(all_candidate_pcs) || ncol(all_candidate_pcs) == 0) {
      warning("No PCs extracted from noise pool. Returning K_star=0.")
      k_actual_max <- 0
      all_candidate_pcs <- NULL # Ensure it's consistently NULL if empty
  } else {
      k_actual_max <- ncol(all_candidate_pcs)
      if (getOption("fmridenoise.verbose_gdlite", FALSE)) message(sprintf("  Extracted %d candidate PCs (up to k_max=%d).", k_actual_max, k_max))
  }
  

  # --- Step 6: Define Good Voxel Mask for K selection --- 
  if (getOption("fmridenoise.verbose_gdlite", FALSE)) message("Step 6: Defining good voxel mask for K selection...")
  if (is_single_run_or_cv_failed) {
    if (getOption("fmridenoise.verbose_gdlite", FALSE)) message("  Single run or LORO CV failed, using all voxels (if any) or TSNR > 0 for good voxel K selection mask.")
    # Fallback for good voxels: use those with TSNR > 0 or simply all voxels if TSNR is also problematic
    good_voxel_mask_for_k_selection <- tsnr_values > 0 
    if (sum(good_voxel_mask_for_k_selection) == 0 && n_voxels > 0) { # If TSNR > 0 yields nothing, take all non-NA TSNR voxels
        good_voxel_mask_for_k_selection <- !is.na(tsnr_values) & tsnr_values != 0 # Avoid all-zero tsnr if possible
         if (sum(good_voxel_mask_for_k_selection) == 0 && n_voxels > 0) good_voxel_mask_for_k_selection <- rep(TRUE, n_voxels) # Absolute fallback
    }
  } else {
    good_voxel_mask_for_k_selection <- (r2_cv_initial >= r2_thresh_good_voxels)
  }
  good_voxel_mask_for_k_selection[is.na(good_voxel_mask_for_k_selection)] <- FALSE # Ensure no NAs
  if (getOption("fmridenoise.verbose_gdlite", FALSE)) message(sprintf("  Good voxel mask for K selection defined with %d voxels.", sum(good_voxel_mask_for_k_selection)))

  # --- Step 7: Select Optimal K_star --- 
  optimal_k_info <- NULL
  if (!is.null(all_candidate_pcs) && k_actual_max > 0) {
      if (getOption("fmridenoise.verbose_gdlite", FALSE)) message("Step 7: Selecting optimal K_star...")
      optimal_k_info <- select_optimal_k_gdlite(
        Y_fmri = Y_fmri,
        X_base_design = X_gd,
        all_pcs = all_candidate_pcs,
        run_idx = run_idx,
        good_voxel_mask = good_voxel_mask_for_k_selection,
        K_max = k_actual_max, # Use actual number of PCs extracted
        min_K = min_K_optimal_selection
      )
  } else {
       if (getOption("fmridenoise.verbose_gdlite", FALSE)) message("Step 7: Skipping optimal K selection as no candidate PCs were extracted.")
      # Create a dummy optimal_k_info if no PCs to test
      optimal_k_info <- list(K_star = 0, selected_pcs = NULL, median_r2cv_values = numeric(0), all_r2cv_per_K = list())
  }
  
  K_star <- optimal_k_info$K_star
  selected_pcs <- optimal_k_info$selected_pcs # This will be NULL if K_star is 0

  # --- Step 8: Final GLM (Optional) --- 
  betas_gdlite <- NULL
  residuals_gdlite <- NULL
  if (perform_final_glm) {
    if (getOption("fmridenoise.verbose_gdlite", FALSE)) message(sprintf("Step 8: Performing final GLM with K_star = %d PCs...", K_star))
    X_final_gdlite <- X_gd
    if (K_star > 0 && !is.null(selected_pcs)) {
      X_final_gdlite <- cbind(X_gd, selected_pcs)
      if (any(duplicated(colnames(X_final_gdlite)))) {
         colnames(X_final_gdlite) <- make.names(colnames(X_final_gdlite), unique=TRUE)
      }
    }
    
    if (ncol(X_final_gdlite) > 0) {
        final_fit <- tryCatch({
            stats::lm.fit(x = X_final_gdlite, y = Y_fmri)
        }, error = function(e) {
            warning(paste("Error in final GLM fit for GLMdenoise-Lite:", e$message))
            return(NULL)
        })
        if(!is.null(final_fit)){
            betas_gdlite <- final_fit$coefficients
            residuals_gdlite <- final_fit$residuals
        }
    } else {
        warning("Final GLMdenoise design matrix X_final_gdlite is empty. Cannot perform final fit.")
    }
  }

  if (getOption("fmridenoise.verbose_gdlite", FALSE)) message("GLMdenoise-Lite Analysis Finished.")

  return(list(
    selected_pcs = selected_pcs,
    K_star = K_star,
    X_gd = X_gd,
    r2_cv_initial = r2_cv_initial,
    tsnr_values = tsnr_values,
    noise_pool_mask = noise_pool_mask,
    good_voxel_mask_for_k_selection = good_voxel_mask_for_k_selection,
    median_r2cv_by_K = optimal_k_info$median_r2cv_values,
    all_r2cv_per_K = optimal_k_info$all_r2cv_per_K,
    betas_gdlite = betas_gdlite,
    residuals_gdlite = residuals_gdlite,
    all_candidate_pcs = all_candidate_pcs # Return all PCs for potential diagnostics
  ))
} 