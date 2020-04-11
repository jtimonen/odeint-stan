#include functions.stan
data {
  int K;                   // number of age classes
  vector[K] age_dist;
  int pop_t;               // total population
  real tswitch;
  
  // Data to fit
  int D;                   // number of days with reported incidence
  int incidence_cases[D];  // overal incidence for W weeks
  int incidence_deaths[D]; // overal incidence for W weeks
  int agedistr_cases[K];   // number of cases at tmax for the K age classes
  int agedistr_deaths[K];  // mortality at tmax for the K age classes
  
  // Prior hyperparameters
  real p_beta;
  real p_eta[2];
  real p_pi[2];
  real p_epsilon[2];
  real p_rho[2];
  real p_phi;
  real p_xi;
  real p_nu;
  real p_psi[2];
  
  // Fixed parameters
  real p_incubation;
  real p_infectious;
  int G;
  real p_gamma[G];
  
  // Simulation
  real t0;           // starting time
  int t_data;        // time of first data
  int S;
  real ts[S];        // time bins
  real contact[K*K];

}

transformed data {
  real EPS = 1.0E-9; // correction to avoid neg. values as tolerance is 1.0E-10
  real x_r[3+K*K+K]; // 4 parameters + K*K contact matrix parameters + K age_dist parameters
  int x_i[1] = {K};
  real init[K*4] = rep_array(0.0, K*4); // initial values
  x_r[1] = tswitch;
  x_r[2] = p_incubation;
  x_r[3] = p_infectious;
  x_r[4:(3+K*K)] = contact;
  for(k in 1:K) {
    x_r[3 + K*K + k] = age_dist[k];
  }
}

parameters{
  real<lower=0,upper=1> beta;             // base transmission rate
  real<lower=0,upper=1> eta;              // reduct. in transm. rate after quarant. measures
  vector<lower=0,upper=1> [K] epsilon;    // age-dependent mortality probability
  vector<lower=0,upper=1> [K-1] raw_rho;  // age-dependent reporting probability
  real<lower=0, upper=1> pii;             // number of cases at t0
  real<lower=0> phi[2];                   // variance parameters
  real<lower=0,upper=1> xi_raw;           // slope of quarantine implementation
  real<lower=0> nu;                       // shift of quarantine implementation
  real<lower=0,upper=1> psi;              // proportion of symptomatics
}

transformed parameters {
  
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
  
  // run ODE solver
  y = integrate_ode_bdf(SEIR, init, t0, ts, theta, x_r, x_i, 1.0E-10, 1.0E-10, 1.0E3);
  
  // Prior
  log_prior_na += log_prior_noadjustment(beta, eta, epsilon, raw_rho, xi_raw, pii, psi, phi, nu,
      p_beta, p_eta, p_epsilon, p_rho, p_xi, p_pi, p_psi, p_phi, p_nu);
  
  // Likelihood
  log_lik_na += log_likelihood_noadjustment(psi, pii, phi, p_gamma, epsilon, rho, incidence_cases, incidence_deaths, agedistr_cases, agedistr_deaths, age_dist, EPS, pop_t, t_data, y);
  
}

model {
  
  // evaluate lp__ with the (invisible) jacobian adjustment term included
  target += log_prior_na;
  target += log_lik_na;
  
}
