% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/900_BTRR_sampling_fns.R
\name{BTRR_draw_sigma_epsilon_sq}
\alias{BTRR_draw_sigma_epsilon_sq}
\title{Draw sigma_epsilon_sq}
\usage{
BTRR_draw_sigma_epsilon_sq(a_epsilon, b_epsilon, Y, x, betas, k)
}
\arguments{
\item{a_epsilon, b_epsilon}{Scalar values of the hyperparameters}

\item{Y}{Array of the response values}

\item{x}{Matrix of the covariates}

\item{betas}{List of length \eqn{D} comprised of \eqn{p_j} by \eqn{R}
matrices}

\item{k}{Scalar for the AR(1) parameter}
}
\value{
A scalar
}
\description{
Draw sigma_epsilon_sq from its posterior full conditional distribution.
}
\keyword{internal}
