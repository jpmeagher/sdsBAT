#' Log Spectral Density Estimation
#'
#' Estimates log10 of the spectral density of a signal.
#'
#' @param signal A univariate time series
#'
#' @return log10 of the spectral density
#'
#' @export
log_spectral_density <- function(signal){
  N <- length(signal)
  spectrum <- fft(signal)
  magnitude_spectrum <- abs(spectrum[1:(floor(N/2)+1)])^2
  log_density <- log10(magnitude_spectrum)

  return(log_density)
}

#' Log Spectral Density Estimation
#'
#' Applies the 'log_spectral_density()' function to the echolocation call
#' dataset
#'
#' @param signals Defaults to the echolocation call dataset provided
#'
#' @return log10 of the spectral density for each call.
#'
#' @seealso \code{\link{log_spectral_density}}
#'
#' @export
log_spectral_density_calls <- function(signals = calls){
  applied_function <- function(signal){
    if (length(signal[[1]]) == 1){
      return(signal[[1]])
    }else{
      return(log_spectral_density(signal[[1]]))
    }
  }
  densities <- apply(signals, c(1,2), applied_function)

  return(as.data.frame(densities))
}

#' Curve Smoothing and Regularisation
#'
#' Fits smoothing splines to noisy observations of a curve and then regularises
#' the curve by interpolating curve values over \eqn{n_points} equidistant
#' points on along the x-axis.
#'
#' @inheritParams log_spectral_density
#' @param observed_range The range over which the noisy input curve is observed at equidistant points.
#' @param regularised_range The range over which the smoothed representation of the curve is to be regularised to.
#' @param n_points Number of equidistant points along the x axis to evaluate the
#'   smoothed curve at.
#'
#' @return A smoothed and regularised representation of the input.
#'
#' @export
smooth_and_regularise <- function(spectrum, observed_range = c(0, 250), regularised_range = c(9,212), n_points = 208){
  N <- length(spectrum)
  observed_frequencies <- seq(observed_range[1], observed_range[2], length.out = N)
  smoothed <- smooth.spline(observed_frequencies, spectrum)
  regularised_frequencies <- seq(regularised_range[1], regularised_range[2], length.out = n_points)
  regularised <- predict(smoothed, regularised_frequencies)$y

  return(regularised)
}

#' Curve Smoothing and Regularisation
#'
#' Applies the 'smooth_and_regularise' function to the spectral densities of the echolocation call
#' dataset
#'
#' @param spectral_densities Defaults to the spectral densities of the echolocation call dataset provided
#'
#' @return smoothed abd regularised spectral density for each call.
#'
#' @seealso \code{\link{log_spectral_density}} \code{\link{log_spectral_density_calls}} \code{\link{smooth_and_regularise}}
#'
#' @export
smooth_and_regularise_call_densities <- function(spectral_densities = log_spectral_density_calls()){
  applied_function <- function(spectral_density){
    if (length(spectral_density[[1]]) == 1){
      return(spectral_density[[1]])
    }else{
      return(smooth_and_regularise(spectral_density[[1]]))
    }
  }
  smooth_regularised <- apply(spectral_densities, c(1, 2), applied_function)

  return(as.data.frame(smooth_regularised))
}

#' Scale curvea
#'
#' Rescales a curve such that it lies along the interval [0, 1].
#'
#' @param spectral_density curve to be mapped to interval
#'
#' @return Curve on interval [0, 1]
#'
#' @export
rescale_density <- function(spectral_density){
  minimum_value <- min(spectral_density)
  scaled_curve <- spectral_density + minimum_value
  new_maximum <- max(scaled_curve)
  scaled_curve <- scaled_curve / new_maximum

  return(scaled_curve)
}

#' Curve Scaling
#'
#' Applies the 'rescale_dinsity' function to the smoothed and regularised
#' spectral densities of the echolocation call dataset.
#'
#' @param spectral_densities Defaults to the smoothed and reguarised spectral
#'   densities of the echolocation call dataset provided.
#'
#' @return smoothed and regularised spectral density for each call scaled to the
#'   interval [0, 1].
#'
#' @seealso \code{\link{log_spectral_density}}
#'   \code{\link{log_spectral_density_calls}}
#'   \code{\link{smooth_and_regularise}}
#'   \code{\link{smooth_and_regularise_densities}} \code{\link{rescale_density}}
#'
#' @export
rescale_smooth_densities <- function(spectral_densities = smooth_and_regularise_call_densities()){
  applied_function <- function(spectral_density){
    if (length(spectral_density[[1]]) == 1){
      return(spectral_density[[1]])
    }else{
      return(rescale_density(spectral_density[[1]]))
    }
  }
  rescaled_densities <- apply(spectral_densities, c(1, 2), applied_function)

  return(as.data.frame(rescaled_densities))
}
