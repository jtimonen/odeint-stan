functions {
#include stan_functions/switch.stan
#include stan_functions/SEIR.stan
#include stan_functions/extract.stan
#include stan_functions/likelihood.stan
#include stan_functions/prior.stan

  // Integrate from one time to another with max timestep dt
  real[] explicit_midpoint_single_output_integrator(real[] init, real t0, real t1, real[] theta, real dt, data real[] x_r, data int[] x_i) {
    vector[size(init)] ytmp = to_vector(init);
    vector[size(init)] midpoint;
    real t = t0;
    while((t + dt) < t1) {
      midpoint = ytmp + (dt / 2.0) * to_vector(SEIR(t, to_array_1d(ytmp), theta, x_r, x_i));
      ytmp = ytmp + dt * to_vector(SEIR(t + dt / 2.0, to_array_1d(midpoint), theta, x_r, x_i));
      t = t + dt;
    }

    midpoint = ytmp + ((t1 - t) / 2.0) * to_vector(SEIR(t, to_array_1d(ytmp), theta, x_r, x_i));
    ytmp = ytmp + (t1 - t) * to_vector(SEIR(t + ((t1 - t) / 2.0), to_array_1d(midpoint), theta, x_r, x_i));
    return to_array_1d(ytmp);
  }

  // Replacement for ode_integrate_rk45 (or ode_integrate_bdf)
  real[,] explicit_midpoint_integrator(real[] init, real t0, real[] ts, real[] theta, real dt, data real[] x_r, data int[] x_i) {
    real y[size(ts), size(init)];
    y[1] = explicit_midpoint_single_output_integrator(init, t0, ts[1], theta, dt, x_r, x_i);
    for(i in 2:size(ts)) {
      y[i] = explicit_midpoint_single_output_integrator(y[i - 1], ts[i - 1], ts[i], theta, dt, x_r, x_i);
    }
    return(y);
  }
}

data {
#include stan_chunks/data.stan
real<lower=0> step_size;
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
  y = explicit_midpoint_integrator(init, t0, ts, theta, step_size, x_r, x_i);
  
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
