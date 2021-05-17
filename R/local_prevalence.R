#' Estimate local prevalence using both targeted and randomised testing data
#'
#' @param test_df A data.frame with the following columns
#'   \itemize{
#'     \item{Nt: }{Total number of targeted tests}
#'     \item{nt: }{Number of targeted tests returning positive}
#'     \item{Nr: }{Total number of randomised tests}
#'     \item{nr: }{Number of randomised tests returning positive}
#'     \item{M: }{Population of local area}
#'     \item{delta_prior_mean: }{Prior mean of bias parameter}
#'     \item{delta_prior_sd: }{Prior standard deviation of bias parameter}
#'     \item{alpha: }{Alpha parameter of beta-binomial prior on probability of infectious given PCR positive}
#'     \item{beta: }{Beta parameter of beta-binomial prior on probability of infectious given PCR positive}
#'   }
#' @param control A list of control parameters defined via
#'   \code{\link{get_control_parameters}}
#' @param imperfect Logical, whether to account for imperfect testing
#' @param type One of PCR_positive or Infectious
#'
#' @export

local_prevalence <- function(test_df,
                             control = get_control_parameters(),
                             imperfect = TRUE,
                             type = "PCR_positive") {

  if (!type %in% c("PCR_positive", "Infectious")) {
    stop("Type must be one of PCR_positive or Infectious")
  }

  # Preallocate log likelihood
  n_bins  <- length(control$I_seq)
  n_time  <- nrow(test_df)
  n_quant <- control$n_quant_approx_bias

  ll_targeted <- array(NA, dim = c(n_time, n_bins, n_quant, 1))

  # Calculate log likelihood
  for (qnum in 1:control$n_quant_approx_bias) {

    test_df$delta <- stats::qnorm(control$quant_approx[qnum],
                                  mean = test_df$delta_prior_mean,
                                  sd = test_df$delta_prior_sd)
    #######################################
    # NOTE: nu currently being fixed -- this may change in future
    #######################################
    test_df$nu <- boot::logit((test_df$Nt - test_df$nt) / test_df$M)

    for (i in seq_along(control$I_seq)) {
      m <- control$I_seq[i]
      which_m_ok <- which(m < test_df$M - (test_df$Nt - test_df$nt))

      if (imperfect) {
        if (length(which_m_ok) > 0) {
          test_df$p1[which_m_ok] <- boot::inv.logit(test_df$nu[which_m_ok] + test_df$delta[which_m_ok])
          test_df$p2[which_m_ok] <- boot::inv.logit(test_df$nu[which_m_ok])

          for (j in seq_along(which_m_ok)) {
            ind <- which_m_ok[j]
            z_eval_points <- stats::qbinom(control$quant_approx, m, test_df$p1[ind])
            z_eval_points <- z_eval_points[z_eval_points <= pmin(m, test_df$Nt[ind])]

            if (length(z_eval_points) > 0) {
              z_eval <- matrix(data = z_eval_points,
                               nrow = length(z_eval_points),
                               ncol = control$n_quant_approx_bias)
              n_type2_eval <- outer(z_eval_points, control$quant_approx,
                                    function(X, Y) stats::qbinom(Y, X, control$beta_testing))
              which_n_type2_eval_ok <- which(n_type2_eval >= pmax(0, z_eval - test_df$nt[ind]) &
                                               n_type2_eval <= pmin(z_eval, test_df$Nt[ind] - test_df$nt[ind]))

              pop_sampling_term_2 <- stats::dbinom(x = test_df$Nt[ind] - c(z_eval),
                                                   size = test_df$M[ind] - m,
                                                   prob = test_df$p2[ind], log = TRUE)
              test_error_term_2 <- stats::dbinom(x = c(n_type2_eval) + test_df$nt[ind] - c(z_eval),
                                                 size = test_df$Nt[ind] - c(z_eval),
                                                 prob = control$alpha_testing, log = TRUE)
              log_summand <- pop_sampling_term_2 + test_error_term_2
              ll_targeted[ind, i, qnum, 1] <- matrixStats::logSumExp(log_summand[which_n_type2_eval_ok])
            }
          }
        }
      } else {
        which_m_ok <- which(m < test_df$M - (test_df$Nt - test_df$nt))

        test_df$p1[which_m_ok] <- boot::inv.logit(test_df$nu[which_m_ok] + test_df$delta[which_m_ok])
        test_df$p2[which_m_ok] <- boot::inv.logit(test_df$nu[which_m_ok])

        ll_targeted_pos <- stats::dbinom(test_df$nt[which_m_ok], m,
                                         test_df$p1[which_m_ok], log = TRUE)
        ll_targeted_neg <- stats::dbinom(test_df$Nt[which_m_ok] - test_df$nt[which_m_ok],
                                         test_df$M[which_m_ok] - m,
                                         test_df$p2[which_m_ok], log = TRUE)
        ll_targeted[which_m_ok, i, qnum, 1] <- ll_targeted_pos + ll_targeted_neg
      }
    }
  }
  ll_targeted[is.na(ll_targeted)] <- -Inf

  # Calculate posterior
  ll_targeted_norm <- apply(ll_targeted, c(1, 3, 4), normalise_logprob)
  ll_prev <- apply(ll_targeted_norm, c(2, 1),
                   function(x) matrixStats::logSumExp(x) - log(n_quant))

  I_prior <- prior_prevalence(test_df, control)
  I_log_post <- ll_prev + log(I_prior)
  I_log_post_max <- apply(I_log_post, 1, max)
  I_post_unnorm <- exp(I_log_post - I_log_post_max)
  I_post_norm <- I_post_unnorm / rowSums(I_post_unnorm)

  I2_post_norm <- I_post_norm
  if (type == "Infectious") {
    # PCR positive -> infectious
    I2_post_norm[] <- 0
    for (i in seq_along(control$I_seq)) {
      for (j in seq_len(nrow(test_df))) {
        bin_probs <- diff(c(0, pmin(cumsum(extraDistr::pbbinom(control$I_seq,
                                                               control$I_seq[i],
                                                               test_df$alpha[j],
                                                               test_df$beta[j])), 1)))
        I2_post_norm[j,] <- I2_post_norm[j,] + (I_post_norm[j,i] * bin_probs)
      }
    }
  }

  I_quant <- t(apply(I2_post_norm, 1, function(v) {
    if(any(is.na(v))) {
      return( rep(NA, control$n_quant_approx_bias))
    } else {
      return(control$I_seq[findInterval(control$quant_approx, cumsum(v)) + 1])
    }
  }))

  list(log_post = I_log_post,
       norm_post = I_post_norm,
       I_quant = I_quant)
}


#' Calculate prior prevalence
#'
#' @inheritParams local_prevalence

prior_prevalence <- function (test_df, control) {

  n_bins  <- length(control$I_seq)
  n_time  <- nrow(test_df)
  I_prior <- matrix(control$bin.d$bin_width, n_time, n_bins,
                    dimnames = list(seq(n_time), control$I_seq),
                    byrow = TRUE)

  for (i in 1:n_time) {
    I_prior[i, ] <- diff(c(0, stats::pnorm(control$bin.d$e.u / test_df$M[i],
                                           mean = control$priors$trunc_gauss$mean,
                                           sd = control$priors$trunc_gauss$sd_proportion)))[1:n_bins]
    I_prior[i, control$I_seq / test_df$M[i] > control$priors$trunc_gauss$upper_trunc_proportion] <- 0
    I_prior[i, "0"] <- I_prior[i, "1"]
  }

  I_prior <- I_prior / rowSums(I_prior)
  I_prior
}


#' Calculate log posterior of prevalence based on data from randomised testing
#'
#' @param test_df A data.frame with the following columns
#'   \itemize{
#'     \item{Nr: }{Total number of randomised tests}
#'     \item{nr: }{Number of randomised tests returning positive}
#'     \item{M: }{Population of local area}
#'   }
#' @param control A list of control parameters defined via
#'   \code{\link{get_control_parameters}}
#' @param imperfect Logical, whether to account for imperfect testing

randomised_testing_prevalence <- function(test_df, control, imperfect) {

  # Preallocate log likelihood
  n_bins  <- length(control$I_seq)
  n_time  <- nrow(test_df)
  ll_randomised <- matrix(NA, n_time, n_bins)

  # Calculate log likelihood
  for (i in seq_along(control$I_seq)) {
    I_curr <- control$I_seq[i]

    if (imperfect) {
      which_I_ok <- which(I_curr < test_df$M)

      if (length(which_I_ok) > 0) {
        for (j in seq_along(which_I_ok)) {

          ind <- which_I_ok[j]
          z_eval_points <- stats::qhyper(p = control$quant_approx, I_curr,
                                         test_df$M[ind] - I_curr,
                                         test_df$Nr[ind])
          z_eval_points <- z_eval_points[z_eval_points <= pmin(I_curr, test_df$Nr[ind])]

          if (length(z_eval_points) > 0) {

            z_eval <- matrix(data = z_eval_points, nrow = length(z_eval_points),
                             ncol = control$n_quant_approx_bias)
            n_type2_eval <- outer(z_eval_points, control$quant_approx,
                                  function(X, Y) stats::qbinom(p = Y, size = X,
                                                               prob = control$beta_testing))

            which_n_type2_eval_ok <- which(n_type2_eval >= pmax(0, z_eval - test_df$nr[ind]) &
                                             n_type2_eval <= pmin(z_eval, test_df$Nr[ind] - test_df$nr[ind]))

            test_error_term2 <- stats::dbinom(c(n_type2_eval) + test_df$nr[ind] - c(z_eval),
                                               test_df$Nr[ind] - c(z_eval),
                                               control$alpha_testing, log = TRUE)

            ll_randomised[ind, i] <- matrixStats::logSumExp(test_error_term2[which_n_type2_eval_ok])
          }
        }
      }
    } else {
      ll_randomised[ ,i] <- stats::dhyper(test_df$nr, I_curr,
                                          pmax(0, test_df$M - I_curr),
                                          test_df$Nr, log = TRUE)
    }
  }

  # Calculate posterior
  unnorm_lik <- exp(ll_randomised - apply(ll_randomised, 1, max))
  norm_lik <- unnorm_lik / rowSums(unnorm_lik)

  I_prior <- prior_prevalence(test_df, control)
  I_log_post <- ll_randomised + log(I_prior)
  I_log_post_max <- apply(I_log_post, 1, max)
  I_post_unnorm <- exp(I_log_post - I_log_post_max)
  I_post_norm <- I_post_unnorm / rowSums(I_post_unnorm)

  post_quant <- t(apply(I_post_norm, 1, function(v) {
    if(any(is.na(v))) {
      return( rep(NA, control$n_quant_approx_bias))
    } else {
      return(control$I_seq[findInterval(control$quant_approx, cumsum(v)) + 1])
    }
  }))

  post_quant
}


#' Calculate log likelihood for targeted testing data
#'
#' @param test_df A data.frame with the following columns
#'   \itemize{
#'     \item{Nt: }{Total number of targeted tests}
#'     \item{nt: }{Number of targeted tests returning positive}
#'     \item{M: }{Population of local area}
#'     \item{p1: }{Probability of taking a test given infectious}
#'     \item{p2: }{Probability of taking a test given not infectious}
#'   }
#' @param I Prevalence
#' @param control A list of control parameters defined via
#'   \code{\link{get_control_parameters}}
#' @param imperfect Logical, whether to account for imperfect testing
#'
#' @export

targeted_testing_loglik <- function(test_df, I, control, imperfect) {

  # Preallocate log likelihood
  n_time  <- nrow(test_df)
  ll_targeted <- double(n_time)

  # Calculate log likelihood
  which_m_ok <- which(I < test_df$M - (test_df$Nt - test_df$nt))

  if (imperfect) {
    if (length(which_m_ok) > 0) {
      for (j in seq_along(which_m_ok)) {
        ind <- which_m_ok[j]
        z_eval_points <- stats::qbinom(control$quant_approx, I[ind], test_df$p1[ind])
        z_eval_points <- z_eval_points[z_eval_points <= pmin(I[ind], test_df$Nt[ind])]

        if (length(z_eval_points) > 0) {
          z_eval <- matrix(data = z_eval_points,
                           nrow = length(z_eval_points),
                           ncol = control$n_quant_approx_bias)
          n_type2_eval <- outer(z_eval_points, control$quant_approx,
                                function(X, Y) stats::qbinom(Y, X, control$beta_testing))
          which_n_type2_eval_ok <- which(n_type2_eval >= pmax(0, z_eval - test_df$nt[ind]) &
                                           n_type2_eval <= pmin(z_eval, test_df$Nt[ind] - test_df$nt[ind]))

          pop_sampling_term_2 <- stats::dbinom(x = test_df$Nt[ind] - c(z_eval),
                                               size = test_df$M[ind] - I[ind],
                                               prob = test_df$p2[ind], log = TRUE)
          test_error_term_2 <- stats::dbinom(x = c(n_type2_eval) + test_df$nt[ind] - c(z_eval),
                                             size = test_df$Nt[ind] - c(z_eval),
                                             prob = control$alpha_testing, log = TRUE)
          log_summand <- pop_sampling_term_2 + test_error_term_2
          ll_targeted[ind] <- matrixStats::logSumExp(log_summand[which_n_type2_eval_ok])
        }
      }
    }
  } else {
    ll_targeted_pos <- stats::dbinom(test_df$nt[which_m_ok], I[which_m_ok],
                                     test_df$p1[which_m_ok], log = TRUE)
    ll_targeted_neg <- stats::dbinom(test_df$Nt[which_m_ok] - test_df$nt[which_m_ok],
                                     test_df$M[which_m_ok] - I[which_m_ok],
                                     test_df$p2[which_m_ok], log = TRUE)
    ll_targeted[which_m_ok] <- ll_targeted_pos + ll_targeted_neg
  }

  ll_targeted[is.na(ll_targeted)] <- -Inf

  ll_targeted
}
