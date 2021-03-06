% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stan_bsrmmlogn.R
\name{stan_bsrmmlogn}
\alias{stan_bsrmmlogn}
\title{Bayesian Structural Risk Measurement Model - Log Normal likelihood}
\usage{
stan_bsrmmlogn(
  losses,
  sbu,
  time,
  bu = NULL,
  n_t = NULL,
  n_sbu = NULL,
  n_bu = NULL,
  frequency_prior = list(a_lambda, b_lambda, r_diric),
  severity_prior = list(scale_sigma = 10, scale_intercept = 10, m0 =
    length(unique(sbu)), sigma_tilde = NULL, slab_scale = 25, slab_df = 25),
  H = 1,
  pars = c("lambdas", "mu", "sigma", "Y_hat"),
  ...
)
}
\arguments{
\item{losses}{a sequence of losses}

\item{sbu}{a sequence of integers indicating which sbu generates a loss.  This
should be of the same length as \code{losses}}

\item{bu}{an optional sequence of the same length as n_sbu}

\item{n_t}{an optional integer allows for manually setting the number of time periods
in case the last time period has no loss in the data.}

\item{n_sbu}{an optional integer allows for manually setting number of sbu in case
the last sbu has no loss in the data.}

\item{n_bu}{an optional integer allows for manually setting number of bu in case
the last bu has no loss in the data.}

\item{frequency_prior}{a list of three prior parameters that characterize
frequency prior. \code{a_lambda} and \code{b_lambda} are unit-length numerics
determining the Gamma prior on the expected frequency at the enterprise level.
\code{r_diric} is a sequence of the same length as the number of sbu, determining
how the frequency at the enterprise level is distributed across sbu.}

\item{severity_prior}{a list of prior parameters for the prior distributions
on the parameters of the log normal likelihood.  \code{a_nu} and \code{b_nu}
are numerics of length 1.

\code{scale_sigma} is the prior parameter on the scale parameter of
half-Cauchy prior on the scale parameter of the log normal likelihood.

\deqn{sigma ~ halfCauchy(scale_sigma)}

Default to \code{scale_sigma} = 10}

\item{pars}{an optional string sequence in \code{rstan::sampling} to select
parameters to extract.  \code{lambdas}, \code{mu} and \code{sigma}
are required to construct the posterior distribution as well as the
simulated loss distribution based on the posterior distribution of the
parameters.}
}
\value{
a \code{stanfit} object output from \code{rstan::sampling}
}
\description{
Fit Bayesian Structural Risk Measurement Model using stan.
This function is a wrapper function around \emph{bstrmmlogn.stan}.
}
\examples{
\dontrun{
# simulate data
nPeriods <- 20
riskProf <- list(
sbu1 = list(lambda = 5, mu = 12, sigma = 2, nu = 10, sbu = 1, bu =1),
sbu2 = list(lambda = 12, mu = 12, sigma = 2, nu = 40, sbu = 2, bu = 1),
sbu3 = list(lambda = 15, mu = 10, sigma = 2, nu = 40, sbu = 3, bu = 1),
sbu4 = list(lambda = 40, mu = 9, sigma = 1.5, nu = 1000, sbu = 4, bu = 2),
sbu5 = list(lambda = 20, mu = 10, sigma = 1.5, nu = 1000, sbu = 5, bu = 2)
)
lossData <- purrr::map_dfr(riskProf, function(l) {
simulateOpLosses(
nPeriods,
lambda = l$lambda,
mu = l$mu,
sigma = l$sigma,
nu = l$nu,
H = 5000
) \%>\% mutate(
sbu = l$sbu,
bu = l$bu
)
})
#fit data with bsrmm with log-normal likelihood
fit <- stan_bsrmmlogn(losses  = lossData$loss_amount,
sbu = lossData$sbu,
time = lossData$time_period,
bu = c(1, 1, 1, 2, 2),
n_t = 20,
H = 5000,
frequency_prior = list(
a_lambda = 10, b_lambda = 4, r_diric = rep(0.2, 5)
),
severity_prior = list(
scale_sigma = 10,
scale_intercept = 10,
m0 = 3,
slab_scale = 10,
slab_df = 20
), sseed = 42, iter = 2000, warmup = 1000, thin = 2, chains = 2
)}

}
