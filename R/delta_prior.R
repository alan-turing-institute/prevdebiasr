#' Specify prior on delta using targeted and randomised regional testing data
#'
#' @param test_df A data.frame with the following columns
#'   \itemize{
#'     \item{Nt: }{Total number of targeted tests in region}
#'     \item{nt: }{Number of targeted tests returning positive in region}
#'     \item{Nr: }{Total number of randomised tests in region}
#'     \item{nr: }{Number of randomised tests returning positive in region}
#'     \item{M: }{Population of region}
#'   }
#' @param control A list of control parameters defined via
#'   \code{\link{get_control_parameters}}
#' @param imperfect Logical, whether to account for imperfect testing
#'
#' @export

specify_delta_prior <- function(test_df,
                                control,
                                imperfect = TRUE) {

  ########################################
  # Calculating I posterior from REACT here
  ########################################
  I_quant <- randomised_testing_prevalence(test_df, control, imperfect)

  ########################################
  # Calculating regional delta posterior
  ########################################
  delta_regional <- delta_regional_posterior(test_df, I_quant, control, imperfect)

  #########################################
  # Specifying smooth delta prior
  #########################################
  n_time <- nrow(test_df)
  ar_cov <- ar1_covariance(
    n_time, control$delta_AR_rho,
    control$delta_AR_sd,
    control$delta_inter_sd
  )

  posterior_mean <- delta_regional$delta_post_mean
  posterior_sd <- delta_regional$delta_post_sd
  posterior_precision <- diag(1 / posterior_sd^2)
  delta_prior_cov <- solve(posterior_precision + solve(ar_cov))
  delta_prior_mean <- c(delta_prior_cov %*% (posterior_precision %*% posterior_mean))
  delta_prior_sd <- sqrt(diag(delta_prior_cov))

  data.frame(
    delta_prior_mean = delta_prior_mean,
    delta_prior_sd = delta_prior_sd
  )
}


#' Approximate delta posterior based on regional testing data
#'
#' @param test_df A data.frame with the following columns
#'   \itemize{
#'     \item{Nt: }{Total number of targeted tests in region}
#'     \item{nt: }{Number of randomised tests returning positive in region}
#'     \item{M: }{Population of region}
#'   }
#' @param control A list of control parameters defined via
#'   \code{\link{get_control_parameters}}
#' @param imperfect Logical, whether to account for imperfect testing
#' @param I_quant Posterior quantiles of prevalence calculated by
#'   \code{\link{randomised_testing_prevalence}}

delta_regional_posterior <- function(test_df, I_quant, control, imperfect) {
  n_time <- nrow(test_df)
  n_quant <- control$n_quant_approx_bias
  n_quant_p1 <- 100

  quant_approx <- seq(n_quant) / (n_quant + 1)
  delquants <- NULL

  # NOTE: fixing nu -- may be changed in future versions of the package
  nu <- boot::logit((test_df$Nt - test_df$nt) / test_df$M)

  test_df$p1_low <- NA
  test_df$p1_high <- NA
  test_df$p2 <- boot::inv.logit(nu)

  for (qnum in 1:n_quant) {
    I <- I_quant[, qnum]

    if (imperfect) {

      # Get range of feasible p1 = ilogit(nu + delta) (uniform prior over this range)
      for (i in 1:n_time) {

        # No false negatives
        true_pos_fn <- min(I[i], stats::qbinom(
          1e-10,
          test_df$nt[i],
          1 - control$alpha_testing
        ))
        test_df$p1_low[i] <- stats::qbeta(
          1e-10,
          true_pos_fn + 1,
          I[i] - true_pos_fn + 1
        )

        # No false positives
        true_pos_fp <- min(I[i], test_df$nt[i] + stats::qbinom(
          1 - 1e-10,
          test_df$Nt[i] - test_df$nt[i],
          control$beta_testing
        ))
        test_df$p1_high[i] <- stats::qbeta(
          1 - 1e-10,
          true_pos_fp + 1,
          I[i] - true_pos_fp + 1
        )
      }

      ll_targeted <- matrix(NA, n_time, n_quant_p1)
      p1_quants <- mapply(
        function(x, y) seq(x, y, length.out = n_quant_p1),
        test_df$p1_low, test_df$p1_high
      )

      for (p1_qnum in 1:n_quant_p1) {
        test_df$p1 <- p1_quants[p1_qnum, ]

        ll_targeted[, p1_qnum] <- targeted_testing_loglik(
          test_df, I,
          control, imperfect
        )
      }
      ll_prev <- apply(ll_targeted, 1, normalise_logprob)

      # Moment matching for beta shape parameters
      p1_mean <- colSums(exp(ll_prev) * p1_quants)
      p1_var <- colSums(exp(ll_prev) * (p1_quants^2)) - (p1_mean^2)
      p1_var <- pmax(1e-10, p1_var) # for numerical stability
      beta_shape1 <- p1_mean * (((p1_mean * (1 - p1_mean)) / p1_var) - 1)
      beta_shape2 <- (1 - p1_mean) * (((p1_mean * (1 - p1_mean)) / p1_var) - 1)

      # Get quantiles for gamma = nu + delta
      gamma_quant <- matrix(NA, n_time, n_quant)
      for (i in which(!is.na(test_df$nt) & (test_df$nt < ((1 - control$beta_testing) * I)))) { # Ignore I quantiles incompatible with nt here
        gamma_quant[i, ] <- boot::logit(stats::qbeta(
          p = quant_approx,
          shape1 = beta_shape1[i],
          shape2 = beta_shape2[i]
        ))
      }
      delta_quant <- gamma_quant - nu
      delquants <- cbind(delquants, delta_quant)
    } else {
      # Below posterior is based on a uniform, i.e. beta(1,1), prior on
      # logit.inv(gamma) = logit.inv(nu + delta)
      beta_shape1 <- test_df$nt + 1
      beta_shape2 <- I - test_df$nt + 1

      gamma_quant <- matrix(NA, n_time, n_quant)
      for (i in which(!is.na(test_df$nt) & beta_shape2 > 0)) { # Ignore I quantiles incompatible with nt here
        gamma_quant[i, ] <- boot::logit(stats::qbeta(
          p = quant_approx,
          shape1 = beta_shape1[i],
          shape2 = beta_shape2[i]
        ))
      }
      delta_quant <- gamma_quant - nu
      delquants <- cbind(delquants, delta_quant)
    }
  }

  delta_post_moment1 <- rowMeans(delquants, na.rm = TRUE)
  delta_post_moment2 <- rowMeans(delquants^2, na.rm = TRUE)
  delta_post_mean <- delta_post_moment1
  delta_post_sd <- sqrt(delta_post_moment2 - (delta_post_moment1)^2)
  
  # When no REACT/ONS, i.e. Nr = 0, output flat prior on delta
  delta_post_mean[test_df$Nr == 0] <- 5
  delta_post_sd[test_df$Nr == 0] <- control$delta_inter_sd
  
  data.frame(delta_post_mean = delta_post_mean, delta_post_sd = delta_post_sd)
}


#' Calculate AR1 covariance matrix
#'
#' @param n_time Number of time points
#' @param rho AR1 correlation term
#' @param noise_sd Standard deviation of noise process
#' @param intercept_sd Standard deviation of intercept

ar1_covariance <- function(n_time, rho, noise_sd, intercept_sd) {

  R_AR1 <- rho ^ abs(outer(1:n_time, 1:n_time, "-")) # AR1 temporal correlation
  R_diagonal <- diag(n_time)
  R_ones <- matrix(1, n_time, n_time)
  delta_cov <- noise_sd^2 * R_AR1 + intercept_sd^2 * R_ones
  return(delta_cov)
}
