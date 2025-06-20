% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ndx_whitening.R
\name{.fit_ar_single_voxel}
\alias{.fit_ar_single_voxel}
\title{Fit AR coefficients for a single voxel}
\usage{
.fit_ar_single_voxel(voxel_residuals, order = 2L, weights_vec = NULL)
}
\arguments{
\item{voxel_residuals}{Numeric vector of residuals for one voxel.}

\item{order}{Integer AR model order.}

\item{weights_vec}{Optional numeric vector of weights for weighted fitting.}
}
\value{
Numeric vector containing the estimated AR coefficients followed by
  the innovation variance.
}
\description{
Helper function used by `ndx_ar2_whitening` to estimate AR coefficients and
innovation variance for one voxel.
}
\keyword{internal}
