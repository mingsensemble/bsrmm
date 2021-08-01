#' Truncated Log Generalized t Random Variables
#'
#' @description Simulate from left-truncated generalized log t distribution
#'
#'@importFrom gamlss.dist qGT
#'@importFrom gamlss.dist pGT
#'@importFrom gamlss.dist rGT
#'
#' @param n an integer giving the number of simulations
#' @param mu a numeric value giving the location parameter for the log
#' generalized t distribution
#' @param sigma a numeric value giving the shape parameter for the log
#' generalized t distribution
#' @param nu a numeric value giving the degree-of-freedom parameter for the log
#' generalized t distribution
#' @param H a numeric value giving the reporting threshold; default to 0.
#' @return a sequence of numerical values of length \code{n} of random draws from
#' Log Generalized t (\code{mu}, \code{sigma}, \code{nu})
#' @examples
#' rLogGtLb(n = 10, mu = 10, sigma = 2, nu = 100, H = 0)
#' rLogGtLb(n = 100, mu = 12, sigma = 2, nu = 20, H = 5000)
#' @export
rLogGtLb <- function(n, mu, sigma, nu, H=0) {
  if(n == 0) {
    draw = 0
  } else {
    if(H > 0) {
      logH <- log(H)
      u <- runif(n)
      x <- gamlss.dist::qGT(
        gamlss.dist::pGT(logH, mu, sigma, nu) +
          u * (1 - gamlss.dist::pGT(logH, mu, sigma, nu)),
        mu, sigma, nu)
    } else {
      # when truncation threshold is 0, sample from full t distribution
      x <- gamlss.dist::rGT(n, mu, sigma, nu)
    }
    draw <- exp(x)
  }
  return(draw)
}

#' Loss Distribution Approach
#'
#' @description This function simulates operational losses following a commonly used loss
#' distribution approach.  Loss counts for each time period is simulated from
#' Poisson distribution with \code{lambda}.  Each loss amount per an loss event
#' is simulated from a log generalized t distribution.
#'
#'@importFrom dplyr bind_cols
#'
#' @param nPeriods An integer giving the number of periods to simulate
#' the losses
#' @param lambda An integer giving the average number of counts per time period
#' @inheritParams rLogGtLb
#' @importFrom gamlss.dist rGT
#' @importFrom dplyr bind_cols
#' @return tibble data frame of loss data at event level.  The first column
#' indicates the time period; the second column indicates the loss amount
#' @examples
#' simulateOpLosses(nPeriods = 10, lambda = 5, mu = 10, sigma = 2, nu = 20, H = 5000)
#' @export
simulateOpLosses <- function(nPeriods, lambda, mu, sigma, nu, H=0, seed = 42) {
  set.seed(seed)

  stopifnot(all(length(lambda) == 1, length(mu) == 1, length(sigma) == 1, length(nu) == 1))

  simCounts <- rpois(nPeriods, lambda)
  simLosses <- do.call('c', sapply(simCounts, function(n){
    rLogGtLb(n, mu, sigma, nu, H)
    }))
  dplyr::bind_cols(
    "time_period" = rep(seq(1, nPeriods), simCounts),
    "loss_amount" = simLosses[simLosses != 0]
    )
}
