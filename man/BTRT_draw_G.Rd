% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/902_BTRTucker_sampling_fns.R
\name{BTRT_draw_G}
\alias{BTRT_draw_G}
\title{Draw core tensor from PFC}
\usage{
BTRT_draw_G(y_til, betas, X, z, V, sig_y2)
}
\arguments{
\item{y_til}{response minus gam_eta}

\item{betas}{tensor decomposition matrices}

\item{X}{tensor covariate}

\item{z}{global variance}

\item{V}{local variance}

\item{sig_y2}{observational variance}
}
\value{
Update to the core tensor
}
\description{
Draw core tensor from PFC
}
\keyword{internal}
