// Likelihood: log normal
// Priors:
// mu: regularized horseshoe
// sigma: hierarchical 

data {
  int<lower = 0> T; // number of periods
  int<lower = 0> N; // number of losses
  int<lower = 0> J; // number of SBUs
  int<lower = 0> K; // number of BUs
  real<lower = 0> H; // reporting threshold
  int<lower = 1, upper = J> SBU[N]; //SBU index
  int<lower = 1, upper = K> BU[J]; // BU index
  vector<lower = H>[N] Y; // loss
  int CNTS[T, J]; // frequency by SBU

  // prior on the frequency
  // lambda prior at enterprise level
  // lambda ~ G(a, b)
  real<lower=0> a_lambda;
  real<lower=0> b_lambda;
  
  // rho_j ~ dirichlet(r)
  vector<lower = 0>[J] r_diric; 
 
  // prior on the severity
  
  // sigma prior
  // sigma ~ C+(scale_sigma)
  real<lower = 0> scale_sigma;

  // prior on intercept
  // b0 ~ N(0, scale_intercept)
  real<lower = 0> scale_intercept;
  
  // horseshoe prior on slope
  // beta ~ N(0, slab_scale^2 * phi_tilde)
  // phi_tilde = sqrt((c^2 * phi^2)/(c^2 + tau^2 * phi^2))
  // phi^2 ~ C+(1)
  // tau ~ C+(tau0)
  real<lower = 0> sigma_tilde; // pseudo standard deviation
  real m0; // expected number of SBU with different mu
  real slab_scale;    // Scale for large slopes
  
  // c^2 ~IG(0.5 * slab_df, 0.5 * slab_df * sigma_tilde^2)
  real slab_df;      // Effective degrees of freedom for large slopes
}

transformed data {
  // horseshoe priors
  real slab_scale2 = square(slab_scale);
  real half_slab_df = 0.5 * slab_df;

  // data transform
  vector<lower = log(H)>[N] logY = log(Y);
  real<lower = 0> logH = log(H); 
  
  // construct a SBU indicator matrix
  matrix[N, J] X = rep_matrix(0, N, J);
  for(n in 1:N) {
    X[n][SBU[n]] = 1;
  }
}

parameters {
  real alpha; // level for mu
  vector[J] beta_tilde; 
  vector<lower=0>[J] phi;
  real<lower=0> c2_tilde;
  real<lower=0> tau_tilde;
  real<lower=0> sigma[J];
  real<lower=0> sigma_bu[K];
  real<lower=0> lambda;
  simplex[J] theta;
}

transformed parameters {
  vector[J] beta;
  vector[J] mu;

  // auxiliary parameters of length N
  vector[N] mu_aux; 
  vector<lower=0>[N] sigma_aux;
  vector<lower=0, upper=1>[N] pp_aux;

  vector<lower=0>[J] sigma_prod_bu; // scale at bu level is assigned to SBU
  vector<lower=0, upper=1>[J] pp; // probability of not being truncated

  // priors on SBU lambdas
  real lambdas[J];
  for(j in 1:J) {
    lambdas[j] = theta[j] * lambda;
  }

  // horseshoe prior
  {
    real tau0 = (m0 / (J - m0)) * (sigma_tilde / sqrt(1.0 * N));
    real tau = tau0 * tau_tilde; // tau ~ cauchy(0, tau0)

    // c2 ~ inv_gamma(half_slab_df, half_slab_df * slab_scale2)
    // Implies that marginally beta ~ student_t(slab_df, 0, slab_scale)
    real c2 = slab_scale2 * c2_tilde;

    vector[J] phi_tilde =
      sqrt( c2 * square(phi) ./ (c2 + square(tau) * square(phi)) );

    // beta ~ normal(0, tau * phi_tilde)
    beta = tau * phi_tilde .* beta_tilde;
    mu_aux = alpha + X * beta;
    // define mu for easy reporting at SBU level
    mu = alpha + beta;
  }

  // hierarchical priors
  // parameters at loss level
  for(i in 1:N) {
    sigma_aux[i] = sigma[SBU[i]];
	pp_aux[i] = pp[SBU[i]];
  }
  // parameters at SBU level
  for(j in 1:J) {
    sigma_prod_bu[j] = sigma_bu[BU[j]];
    pp[j] = 1 - lognormal_cdf(logH, mu[j], sigma[j]);
  }
}
model {
  vector[N] log_lik;
  vector[J] log_lik_cnts;

  // prior on lambda: poisson gamma 
  lambda ~ gamma(a_lambda, b_lambda);
  theta ~ dirichlet(r_diric);
  
  // hyper priors

  // sigma has a common prior at BU level
  sigma_bu ~ cauchy(0, scale_sigma);
  sigma ~ cauchy(0, sigma_prod_bu);

  // mu has a horseshoe prior
  alpha ~ normal(0, scale_intercept);

  beta_tilde ~ normal(0, 1);
  phi ~ cauchy(0, 1);
  tau_tilde ~ cauchy(0, 1);
  c2_tilde ~ inv_gamma(half_slab_df, half_slab_df);

  // likelihood
  // poisson likelihood
  for(j in 1:J) {
    log_lik_cnts[j] = poisson_lpmf(CNTS[, j]|lambdas[j] * pp[j]);
  }
  // log normal likelihood
  log_lik = lognormal_lpdf(logY |mu_aux, sigma_aux) - log(pp_aux);
  target += sum(log_lik) + sum(log_lik_cnts);
}
generated quantities {
  real<lower=0> Y_hat[J];
  int<lower=0> lambda_hat;
  for(j in 1:J) {
    lambda_hat = poisson_rng(lambdas[j]);
	if(lambda_hat == 0) {
	  Y_hat[j] = 0;
	} else {
	  Y_hat[j] = 0;
	  for(k in 1:lambda_hat) {
	    Y_hat[j] += lognormal_rng(mu[j], sigma[j]);
	  }
	}
  }
}
