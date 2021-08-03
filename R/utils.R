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
    out <- values[findInterval(quant_approx, cumsum(x)) + 1]
  }

  out
}


#' Gaussian moment-matching on logit scale
#'
#' @param x Vector of log-probabilities
#' @param values Vector of proportions
#'
#' @export

logit_mom <- function(x, values) {

  if (any(values < 0))
    stop("Proportions must be non-negative")

  logit_values <- boot::logit(values)

  norm_prob <- exp(x) / sum(exp(x))

  # Place mass on zero on smallest non-zero value
  norm_prob_dezeroed <- norm_prob
  norm_prob_dezeroed[1] <- 0
  norm_prob_dezeroed[2] <- norm_prob[1] + norm_prob[2]

  mom_match1 <- sum(norm_prob_dezeroed * logit_values, na.rm = TRUE)
  mom_match2 <- sum(norm_prob_dezeroed * (logit_values ^ 2), na.rm = TRUE)

  data.frame(mean = mom_match1,
             sd =  sqrt(mom_match2 - (mom_match1 ^ 2)))
}

