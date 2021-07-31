#' Capital Allocation
#'
#' @description Calculate capital allocation using Euler Principle
#' based on expected shortfall measures
#' @importFrom rstan extract
#' @param stanfit a \code{stanfit} object of outputs from \code{stan_bsrmmlogt}
#' or \code{stan_bsrmmlogn}
#' @param p percentile for expected shortfall
#' @return a list with two elements:\code{risk measures} capital charges
#' for each segment and \code{allocation} percentage of capital allocation
calculate_allocation <- function(stanfit, ...) {
  UseMethod("calculate_allocation")
}

calculate_allocation.stanfit <- function(stanfit, p) {
  sim_res <- rstan::extract(stanfit, pars = "Y_hat")$Y_hat
  z <- rowSums(sim_res)
  var_hat <- quantile(z, p)
  risk_measures <- colSums(sim_res[(z >= var_hat), ])/sum((z >= var_hat))
  allocated_pp <- risk_measures/sum(risk_measures) * 100
  return(list(
    `risk measures` = risk_measures,
    `allocation` = allocated_pp
  ))
}
