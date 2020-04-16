functions {
#include stan_functions/switch.stan
#include stan_functions/SEIR.stan
#include stan_functions/extract.stan
#include stan_functions/likelihood.stan
#include stan_functions/prior.stan

  // Vector version of SEIR
  vector odefun(real t, vector y, real[] theta, data real[] x_r, data int[] x_i){
    return to_vector(SEIR(t, to_array_1d(y), theta, x_r, x_i));
  }
  
  // Linear interpolation from equispaced time grid
  real[,] interpolate(vector[] y, data int[] R, data real[] A){
    int N = size(R);
    int D = num_elements(y[1]);
    real y_out[N, D];
    for(i in 1:N){
      int idx = R[i];
      y_out[i] = to_array_1d(A[i] * y[idx] + (1 - A[i]) * y[idx+1]);
    }
    return y_out;
  }
  
  // Euler method
  real[,] integrate_ode_euler(real[] y0, real t0, real[] ts, 
      real[] theta, data real h, data int[] INTERP_R, data real[] INTERP_A, 
      data real[] x_r, data int[] x_i){
    int D = size(y0);
    int N = size(ts);
    int R_max = INTERP_R[N];
    real y_out[N, D];
    real t_tmp = t0;
    vector[D] y_tmp = to_vector(y0);
    vector[D] y[R_max+1];
    y[1] = y_tmp;
    for(i in 1:R_max){
      y_tmp = y_tmp + h*odefun(t_tmp, y_tmp, theta, x_r, x_i);
      t_tmp = t_tmp + h;
      y[i+1] = y_tmp;
    }
    y_out = interpolate(y, INTERP_R, INTERP_A);
    return(y_out);
  }
  

}

data {
#include stan_chunks/data.stan
real<lower=0> step_size;
int<lower=0> INTERP_R[D];
real<lower=0> INTERP_A[D];
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
  y = integrate_ode_euler(init, t0, ts, theta, step_size, INTERP_R, INTERP_A, x_r, x_i);
  
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

