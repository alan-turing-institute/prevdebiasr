#' Normalise a vector of log probabilities
#'
#' @param x A vector of log probabilities

normalise_logprob <- function(x) {
  lognorm <- matrixStats::logSumExp(x)
  x - lognorm
}
