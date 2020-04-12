
  // these store non-jacobian-adjusted versions of prior and likelihood
  real log_prior_na = 0.0;
  real log_lik_na = 0.0;
  
  // transformed parameters
  vector[K] rho;
  real xi = xi_raw+0.5;
  real theta[6];          // vector of parameters
  real y[S,K*4];          // raw ODE output
  
  // transformed paremeters
  for(i in 1:(K-1)){
    rho[i] = raw_rho[i];
  }
  rho[K] = 1.0;
  theta[1:6] = {beta, eta, xi, nu, pii, psi};
  