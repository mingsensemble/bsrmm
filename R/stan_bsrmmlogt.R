#' Bayesian Structural Risk Measurement Model - Log t likelihood
#'
#' @description Fit Bayesian Structural Risk Measurement Model using stan.
#' This function is a wrapper function around *bstrmmlogt.stan*.
#'
#' @importFrom loo nlist
#' @import tidyr
#' @import dplyr
#' @param losses a sequence of losses
#' @param sbu a sequence of integers indicating which sbu generates a loss.  This
#' should be of the same length as \code{losses}
#' @param bu an optional sequence of the same length as n_sbu
#' @param n_sbu an optional integer allows for manually setting number of sbu in case
#' the last sbu has no loss in the data.
#' @param n_bu an optional integer allows for manually setting number of bu in case
#' the last bu has no loss in the data.
#' @param n_t an optional integer allows for manually setting the number of time periods
#' in case the last time period has no loss in the data.
#' @param frequency_prior a list of three prior parameters that characterize
#' frequency prior. \code{a_lambda} and \code{b_lambda} are unit-length numerics
#' determining the Gamma prior on the expected frequency at the enterprise level.
#' \code{r_diric} is a sequence of the same length as the number of sbu, determining
#' how the frequency at the enterprise level is distributed across sbu.
#' @param severity_prior a list of prior parameters for the prior distributions
#' on the parameters of the log generalized t likelihood.  \code{a_nu} and \code{b_nu}
#' are numerics of length 1.  \code{a_nu} and \code{b_nu} determines the Gamma
#' prior on the degree-of-freedom parameter \code{nu}.
#'
#' \deqn{nu ~ G(a_nu, b_nu)}
#'
#' Default to a weakly informative prior with \code{a_nu} = 3 and \code{b_nu} = 0.1.
#'
#' \code{scale_sigma} is the prior parameter on the scale parameter of
#' half-Cauchy prior on the scale parameter of the log generalized t likelihood.
#'
#' \deqn{sigma ~ halfCauchy(scale_sigma)}
#'
#' Default to \code{scale_sigma} = 10
#'
#'@param pars an optional string sequence in \code{rstan::sampling} to select
#'parameters to extract.  \code{lambdas}, \code{mu}, \code{sigma}, and \code{nu}
#'are required to construct the posterior distribution as well as the
#'simulated loss distribution based on the posterior distribution of the
#'parameters.
#'
#' @export
stan_bsrmmlogt <- function(losses, sbu, time, bu = NULL, n_t = NULL,
                         n_sbu = NULL, n_bu = NULL,
                         frequency_prior = list(
                           a_lambda, b_lambda, r_diric
                         ),
                         severity_prior = list(
                           a_nu = 3, b_nu = 0.1,
                           scale_sigma = 10,
                           scale_intercept = 10,
                           m0 = length(unique(sbu)),
                           sigma_tilde = NULL,
                           slab_scale = 25, slab_df = 25
                         ), H = 1, pars = c("lambdas", "mu", "sigma", "nu", "Y_hat"),
                         ...
) {

  if(any(!c("lambdas", "mu", "sigma", "nu") %in% pars)) {
    stop("pars have to consists of all of lambdas, mu, sigma, nu to construct posterior distribution")
  }

  n_sbu <- ifelse(is.null(n_sbu), max(sbu), n_sbu)
  n_bu <- ifelse(is.null(n_bu), max(bu), n_bu)
  n_t <- ifelse(is.null(n_t), max(time) - min(time) + 1, n_t)
  if(is.null(bu)) {
    bu <- unique(sbu)
  }

  if(is.null(severity_prior$sigma_tilde)) {
    severity_prior$sigma_tilde <- sd(log(losses))
  }

  # error handling
  stopifnot(length(frequency_prior$r_diric) == n_sbu)
  if(sum(frequency_prior$r_diric) != 1) {
    frequency_prior$r_diric <- frequency_prior$r_diric/sum(frequency_prior$r_diric)
    #r_diric sums up to 1
  }
  stopifnot(length(sbu) == length(losses))
  stopifnot(length(bu) == n_sbu)

  dfBase <- expand.grid(seq(n_t), seq(n_sbu))
  names(dfBase) <- c("time", "sbu")
  lossCounts <- aggregate(losses, list(time, sbu), length)
  names(lossCounts) <- c("time", "sbu", "counts")

  cnts <- dfBase %>% dplyr::left_join(lossCounts, by = c("time", "sbu")) %>%
    tidyr::replace_na(list(counts = 0)) %>% tidyr::pivot_wider(
      id_cols = 'time',
      names_from = "sbu",
      values_from = "counts"
      ) %>% dplyr::select(-time)

  data_input <- list(
    T = n_t,
    N = length(losses),
    J = n_sbu,
    K = n_bu,
    H = H,
    CNTS = cnts
  )

  data_list <- loo::nlist(losses, sbu, bu)
  names(data_list) <- c("Y", "SBU", "BU")

  standata <- c(data_input, data_list, frequency_prior, severity_prior)

  out <- rstan::sampling(stanmodels$bsrmmlogt, data = standata, pars = pars, ...)
  return(out)

}
