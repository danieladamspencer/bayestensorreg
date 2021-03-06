% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/901_BTRR_GGM_sampling_fns.R
\name{BTRR_GGM_extract_last_iter}
\alias{BTRR_GGM_extract_last_iter}
\title{Extract the last iteration of the MCMC from a result file}
\usage{
BTRR_GGM_extract_last_iter(result_file)
}
\arguments{
\item{result_file}{A character string giving the relative file path to a
desired result file.}
}
\value{
A list of the elements from the most recent MCMC draw
}
\description{
Extract the last iteration of the MCMC from a result file
}
\keyword{internal}
