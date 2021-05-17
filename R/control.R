#' Set control parameters
#'
#' @param delta_AR_rho Correlation in AR1 component of the temporal Gaussian process on delta
#' @param delta_AR_sd Standard deviation in AR1 component of the temporal Gaussian process on delta
#' @param delta_inter_sd Prior standard deviation on delta intercept
#' @param alpha_testing False positive rate
#' @param beta_testing False negative rate
#' @param I_prior_mean Mean of truncated Gaussian prior on prevalence
#' @param I_prior_sd Standard deviation of truncated Gaussian prior on prevalence
#' @param I_prior_upper Upper limit of truncated Gaussian prior on prevalence
#' @param I_max Maximum value of prevalence considered in the model
#' @param nbin The number of bins into which prevalence I is divided
#' @param n_quant_approx_bias Number of quantiles to use to approximate delta EB posterior
#'
#' @export

get_control_parameters <- function(delta_AR_rho = 0.99,
                                   delta_AR_sd = 1,
                                   delta_inter_sd = 10,
                                   I_prior_mean = 0.005,
                                   I_prior_sd = 0.025,
                                   I_prior_upper = 0.04,
                                   alpha_testing = 0.001,
                                   beta_testing = 0.05,
                                   I_max = 1e6,
                                   nbin = 200,
                                   n_quant_approx_bias = 25
                                   ) {

  ###################################
  # PRIOR HYPERPARAMETERS
  ####################################
  priors <- list()
  priors$trunc_gauss <- list(mean_proportion = I_prior_mean,
                             sd_proportion = I_prior_sd,
                             upper_trunc_proportion = I_prior_upper)

  ###################################
  # QUANTILE APPROXIMATION
  ####################################
  quant_approx <- (1:n_quant_approx_bias) / (n_quant_approx_bias + 1)


  #########################################################
  # Create bins for # infecteds I
  #########################################################

  ########################################################
  # Code for finding eps - it is determined to a narrow interval
  # once defaults$I_max and nbin are specified
  ########################
  eps <- 0.1
  ok <- F
  while(!ok){
    eps <- eps - .000005
    ev <- c()
    ev[1] <- 0
    ev[2] <- 1
    for(en in 3:(nbin + 1)){
      ev[en] <- ceiling(ev[en - 1] * (1 + eps))[]
    }
    ok <- ev[length(ev)] < I_max
  }
  ev[length(ev)] <- I_max + 1
  bin.d <- data.frame(e.l = ev[1:(length(ev) - 1)], e.u = ev[2:length(ev)] - 1)
  bin.d$mu <- floor((bin.d$e.l + bin.d$e.u) / 2)
  bin.d$bin_width <- bin.d$e.u - bin.d$e.l + 1
  mseq <- bin.d$mu

  I_seq <- bin.d$mu

  list(delta_AR_rho = delta_AR_rho,
       delta_AR_sd = delta_AR_sd,
       delta_inter_sd = delta_inter_sd,
       priors = priors,
       n_quant_approx_bias = n_quant_approx_bias,
       quant_approx = quant_approx,
       alpha_testing = alpha_testing,
       beta_testing = beta_testing,
       bin.d = bin.d,
       I_seq = I_seq,
       nbin = nbin
       )
}
