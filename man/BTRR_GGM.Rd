% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/401_BTRR_GGM.R
\name{BTRR_GGM}
\alias{BTRR_GGM}
\title{Bayesian Tensor Response Regression with Gaussian Graphical Model}
\usage{
BTRR_GGM(
  input,
  n_iter = 100,
  n_burn = 0,
  Rank = 1,
  hyperparameters = NULL,
  save_after = NULL,
  save_llik = TRUE,
  results_file = NULL,
  num_cores = parallel::detectCores() - 2,
  save_dir = "."
)
}
\arguments{
\item{input}{an object of class \code{TRR_GGM_data} or a list with elements
\code{Y} (a list of G arrays with dimensions
\eqn{p_1\times \cdots \times p_D \times T \times n}) and
\code{x} (a \eqn{T \times n} matrix).}

\item{n_iter}{(a scalar) the number of posterior samples desired}

\item{n_burn}{(a scalar) the number of posterior samples to discard as a
burn-in}

\item{Rank}{(a positive integer) the rank for the PARAFAC/CP
tensor decomposition}

\item{hyperparameters}{a list with named numbers containing at least one of
the following: \code{a.tau}, \code{b.tau}, \code{a.lambda},
\code{b.lambda}, \code{a.epsilon}, or \code{b.epsilon} defining the values
of the hyperparameters within the model. If \code{NULL}, then default
values will be used.}

\item{save_after}{(an integer) An .rds file will be saved every
\code{save_after} MCMC iterations with all of the results to that point
(helpful for getting intermediate results in less time)}

\item{save_llik}{(a logical) Should the log-likelihood be calculated and
saved at each iteration? Doing so comes at a cost, but can be used for
model diagnostics. Defaults to \code{TRUE}.}

\item{results_file}{(optional) The relative path to a result file. This is
used to continue an MCMC chain on a set of data for which some results
already exist.}

\item{num_cores}{The number of cores used for running the code in parallel}

\item{save_dir}{(a character) A path to a directory in which the temporary
results will be saved. Defaults to the current working directory.}
}
\value{
A list object with the posterior draws from the MCMC chain.
}
\description{
Performs the MCMC to draw from the posterior distribution for the model
published by Spencer, Guhaniyogi, and Prado (2020).
}
\examples{
\dontrun{
input <- TRR_GGM_simulated_data()
results <- BTRR_GGM(input)
}
}
