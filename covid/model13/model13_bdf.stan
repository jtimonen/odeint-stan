functions {
#include stan_functions/switch.stan
#include stan_functions/SEIR.stan
#include stan_functions/extract.stan
#include stan_functions/likelihood.stan
#include stan_functions/prior.stan
}

data {
#include stan_chunks/data.stan
}

transformed data {
#include stan_chunks/transformed_data.stan
}

parameters{
#include stan_chunks/parameters.stan
}

transformed parameters {
  
  // These store non-jacobian-adjusted versions of prior and likelihood
  real log_prior_na = 0.0;
  real log_lik_na = 0.0;
  
  // Transformed parameters
  vector[K] rho;
  real xi = xi_raw+0.5;
  real theta[6];
  real y[S,K*4];
  for(i in 1:(K-1)){ rho[i] = raw_rho[i]; }
  rho[K] = 1.0;
  theta[1:6] = {beta, eta, xi, nu, pii, psi};
  
  // Solve ODE
  y = integrate_ode_bdf(SEIR, init, t0, ts, theta, x_r, x_i, abs_tol, rel_tol, max_iter);
  
  // Compute prior and likelihood
  log_prior_na += log_prior_noadjustment(beta, eta, epsilon, raw_rho, xi_raw, pii, psi, phi, nu, p_beta, p_eta, p_epsilon, p_rho, p_xi, p_pi, p_psi, p_phi, p_nu);
  log_lik_na += log_likelihood_noadjustment(psi, pii, phi, p_gamma, epsilon, rho, incidence_cases, incidence_deaths, agedistr_cases, agedistr_deaths, age_dist, EPS, pop_t, t_data, y);
}

model {
  // Evaluate lp__ with the (invisible) jacobian adjustment term included
  target += log_prior_na;
  target += log_lik_na;
}

generated quantities{
  real log_prior_na_GEN_ = 0.0;
  real log_lik_na_GEN_ = 0.0;
  
  // Solve ODE
  real y_GEN_[S,K*4] = integrate_ode_bdf(SEIR, init, t0, ts, theta, x_r, x_i, abs_tol, rel_tol, max_iter);
  
  // Compute prior and likelihood
  log_prior_na_GEN_ += log_prior_noadjustment(beta, eta, epsilon, raw_rho, xi_raw, pii, psi, phi, nu, p_beta, p_eta, p_epsilon, p_rho, p_xi, p_pi, p_psi, p_phi, p_nu);
  log_lik_na_GEN_ += log_likelihood_noadjustment(psi, pii, phi, p_gamma, epsilon, rho, incidence_cases, incidence_deaths, agedistr_cases, agedistr_deaths, age_dist, EPS, pop_t, t_data, y_GEN_);
}