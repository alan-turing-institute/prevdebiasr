#' Normalise a vector of log probabilities
#'
#' @param x A vector of log probabilities

normalise_logprob <- function(x) {
  lognorm <- matrixStats::logSumExp(x)
  x - lognorm
}


#' Extract quantiles from a (normalised) probability distribution
#'
#' @param x Vector of (normalised) probabilities
#' @param values Vector of corresponding values
#' @param n_quantiles Number of quantiles
#'
#' @export

get_quantiles <- function(x, values, n_quantiles) {
  quant_approx <- (1:n_quantiles) / (n_quantiles + 1)

  if (any(is.na(x))) {
    out <- rep(NA, n_quantiles)
  } else {
    out <- return(values[findInterval(quant_approx, cumsum(v)) + 1])
  }

  out
}
