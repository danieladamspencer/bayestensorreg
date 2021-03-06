% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/902_BTRTucker_sampling_fns.R
\name{BTRT_draw_Bj}
\alias{BTRT_draw_Bj}
\title{Draw from PFC for betas}
\usage{
BTRT_draw_Bj(y_til, X, betas, tau, W, sig_y2, G, j)
}
\arguments{
\item{y_til}{response minus etagam}

\item{X}{covariate tensor}

\item{betas}{all betas, for betas from other dimensions}

\item{tau}{global variance}

\item{W}{local variance}

\item{sig_y2}{observational variance}

\item{G}{core tensor}

\item{j}{dimension index}
}
\value{
updates the betas object at index j
}
\description{
Draw from PFC for betas
}
\keyword{internal}
