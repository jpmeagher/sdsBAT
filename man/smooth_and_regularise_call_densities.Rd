% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/call_to_spectralCurve.R
\name{smooth_and_regularise_call_densities}
\alias{smooth_and_regularise_call_densities}
\title{Curve Smoothing and Regularisation}
\usage{

  smooth_and_regularise_call_densities(spectral_densities = log_spectral_density_calls())
}
\arguments{
\item{spectral_densities}{Defaults to the spectral densities of the echolocation call dataset provided}
}
\value{
smoothed abd regularised spectral density for each call.
}
\description{
Applies the 'smooth_and_regularise' function to the spectral densities of the echolocation call
dataset
}
\seealso{
\code{\link{log_spectral_density}} \code{\link{log_spectral_density_calls}} \code{\link{smooth_and_regularise}}
}
