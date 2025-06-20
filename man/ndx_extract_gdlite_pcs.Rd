% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ndx_gdlite.R
\name{ndx_extract_gdlite_pcs}
\alias{ndx_extract_gdlite_pcs}
\title{Extract GLMdenoise-Lite Principal Components}
\usage{
ndx_extract_gdlite_pcs(
  Y_fmri,
  X_nuisance_base,
  n_pcs,
  voxel_mask = NULL,
  min_variance_voxels = NULL,
  max_variance_voxels_pca_prop = 0.2,
  center = TRUE,
  scale. = FALSE
)
}
\arguments{
\item{Y_fmri}{Numeric matrix of fMRI data (timepoints x voxels).}

\item{X_nuisance_base}{Numeric matrix representing the design matrix of
effects to regress out before PCA (e.g., task effects from a canonical HRF,
motion parameters, polynomial trends). Must have the same number of rows
as `Y_fmri`.}

\item{n_pcs}{Integer, the number of principal components to extract.}

\item{voxel_mask}{Optional. A logical or numeric vector indicating which voxels
(columns of `Y_fmri`) to include in the PCA. If logical, length must match
`ncol(Y_fmri)`. If numeric, contains column indices. If `NULL` (default),
PCA will be performed on a subset of high-variance voxels (see `min_variance_voxels`, `max_variance_voxels_pca_prop`).}

\item{min_variance_voxels}{Integer, the minimum number of voxels to use for PCA
if `voxel_mask` is `NULL`. Defaults to `max(100, n_pcs * 5)`.}

\item{max_variance_voxels_pca_prop}{Numeric, proportion of total voxels to consider as
candidates for high-variance selection if `voxel_mask` is `NULL`. E.g. 0.2 means top 20% variance voxels.
The actual number used will be `min(max_variance_voxels, number_of_voxels_meeting_prop_threshold)`.
Defaults to `0.2`.}

\item{center}{Logical, whether to center the data before PCA (passed to `stats::prcomp`). Defaults to `TRUE`.}

\item{scale.}{Logical, whether to scale the data before PCA (passed to `stats::prcomp`). Defaults to `FALSE`.}
}
\value{
A matrix containing the extracted temporal principal components
  (timepoints x `n_pcs`). Returns `NULL` if critical errors occur (e.g.,
  `n_pcs` is invalid or no valid voxels for PCA).
}
\description{
Performs a basic GLM to obtain residuals from fMRI data after regressing out
a base set of nuisance regressors (e.g., task, motion, trends). Then, PCA is
applied to these residuals to extract temporal principal components, similar to
a simplified GLMdenoise approach.
}
\examples{
\dontrun{
# Basic Example
n_time <- 200
n_voxels <- 500
n_task_regs <- 2
n_motion_regs <- 6

Y_fmri_example <- matrix(rnorm(n_time * n_voxels), n_time, n_voxels)
X_base_example <- matrix(rnorm(n_time * (n_task_regs + n_motion_regs)), n_time,
                        n_task_regs + n_motion_regs)

gdlite_pcs <- ndx_extract_gdlite_pcs(
  Y_fmri = Y_fmri_example,
  X_nuisance_base = X_base_example,
  n_pcs = 5
)

if (!is.null(gdlite_pcs)) {
  print(paste("Extracted GD-Lite PCs dimensions:", dim(gdlite_pcs)[1], "x", dim(gdlite_pcs)[2]))
}

# With a voxel mask
voxel_mask_example <- sample(c(TRUE, FALSE), n_voxels, replace = TRUE, prob = c(0.3, 0.7))
gdlite_pcs_masked <- ndx_extract_gdlite_pcs(
  Y_fmri = Y_fmri_example,
  X_nuisance_base = X_base_example,
  n_pcs = 3,
  voxel_mask = voxel_mask_example
)
if (!is.null(gdlite_pcs_masked)) {
  print(paste("Extracted masked GD-Lite PCs dimensions:",
               dim(gdlite_pcs_masked)[1], "x", dim(gdlite_pcs_masked)[2]))
}
}
}
