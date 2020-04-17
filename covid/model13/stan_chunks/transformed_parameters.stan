  // These store non-jacobian-adjusted versions of prior and likelihood
  real log_prior_na = 0.0;
  real log_lik_na = 0.0;
  
  // Transformed parameters
  real xi = xi_raw+0.5;
  real theta[6];
  real y[S,K*4];
  theta[1:6] = {beta, eta, xi, nu, pii, psi};
