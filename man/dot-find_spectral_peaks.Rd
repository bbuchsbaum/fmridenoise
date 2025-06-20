% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ndx_spectral.R
\name{.find_spectral_peaks}
\alias{.find_spectral_peaks}
\title{Identify prominent spectral peaks}
\usage{
.find_spectral_peaks(
  mt_res,
  TR,
  nyquist_guard_factor,
  n_sine_candidates,
  verbose = FALSE
)
}
\description{
Searches for peaks below a guarded Nyquist frequency and returns the selected
frequencies in Hz.
}
\keyword{internal}
