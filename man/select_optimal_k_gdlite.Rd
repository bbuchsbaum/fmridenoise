% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ndx_gdlite.R
\name{select_optimal_k_gdlite}
\alias{select_optimal_k_gdlite}
\title{Select Optimal Number of Principal Components for GLMdenoise}
\usage{
select_optimal_k_gdlite(
  Y_fmri,
  X_base_design,
  all_pcs,
  run_idx,
  good_voxel_mask,
  K_max = NULL,
  min_K = 0
)
}
\arguments{
\item{Y_fmri}{Numeric matrix of fMRI data (timepoints x voxels), concatenated across runs.}

\item{X_base_design}{Numeric matrix, the base design matrix (e.g., task regressors,
motion parameters, polynomial trends) to which PCs will be appended.
Must have the same number of rows as `Y_fmri`.}

\item{all_pcs}{Numeric matrix (timepoints x `K_max`), containing all candidate PCs
extracted from noise-pool residuals (e.g., output of `ndx_extract_gdlite_pcs`
called with `n_pcs = K_max`).}

\item{run_idx}{Numeric vector indicating run membership for each timepoint in `Y_fmri`.}

\item{good_voxel_mask}{Logical vector of length `ncol(Y_fmri)`, where TRUE indicates
a "good voxel" over which the median LORO R-squared will be calculated to select `K_star`.
These are typically voxels with a reasonable initial R-squared from a simpler model.}

\item{K_max}{Integer, the maximum number of PCs to consider. Should match `ncol(all_pcs)`.
If `NULL`, it's derived from `ncol(all_pcs)`.}

\item{min_K}{Integer, the minimum number of PCs to consider in the loop. Default is 0 (model with no PCs).}
}
\value{
A list containing:
  - `K_star`: The optimal number of PCs selected.
  - `selected_pcs`: The matrix of selected PCs (`all_pcs[, 1:K_star, drop = FALSE]`).
    If `K_star` is 0, this will be `NULL` or an empty matrix.
  - `median_r2cv_values`: A numeric vector of median R-squared_cv values for each `K` tested.
  - `all_r2cv_per_K`: A list where each element is the full r2_cv vector for all voxels for a given K.
}
\description{
This function determines the optimal number of GLMdenoise principal components (PCs)
to include in a GLM. It iterates from 1 to `K_max` PCs, appends them to a
base design matrix (`X_base_design`), and recalculates leave-one-run-out (LORO)
cross-validated R-squared (`r2_cv`). The optimal number of PCs (`K_star`) is
chosen as the one that maximizes the median `r2_cv` across a predefined set of
"good voxels".
}
