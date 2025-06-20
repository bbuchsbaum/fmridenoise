% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ndx_workflow_utils.R
\name{ndx_run_annihilation_setup}
\alias{ndx_run_annihilation_setup}
\title{Run the optional GLMdenoise-Lite setup for Annihilation mode}
\usage{
ndx_run_annihilation_setup(
  Y_fmri,
  events,
  motion_params,
  run_idx,
  TR,
  opts_annihilation,
  verbose = TRUE
)
}
\arguments{
\item{Y_fmri}{Numeric matrix of fMRI data.}

\item{events}{Data frame of events.}

\item{motion_params}{Motion parameters matrix.}

\item{run_idx}{Run index vector.}

\item{TR}{Repetition time in seconds.}

\item{opts_annihilation}{List of annihilation options.}

\item{verbose}{Logical verbosity flag.}
}
\value{
List with `selected_pcs` and full `gdlite_results` if run.
}
\description{
Executes a light-weight GLMdenoise procedure when Annihilation Mode is enabled. The selected principal components are returned along with the full GLMdenoise-Lite results.
}
