% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ndx_utils.R
\name{\%||\%}
\alias{\%||\%}
\title{Null-Coalescing Operator}
\usage{
a \%||\% b
}
\arguments{
\item{a}{First value.}

\item{b}{Fallback value.}
}
\value{
`a` if not `NULL`, otherwise `b`.
}
\description{
Return `b` when `a` is `NULL`, otherwise `a`.
}
\keyword{internal}
