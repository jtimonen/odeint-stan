
functions {

  real switch_zero(real t, real t1) {
    return(1/(1+exp(5*(t-t1))));
  }
  
  real switch_eta(real t, real t1, real eta, real nu, real xi) {
    return(eta+(1-eta)/(1+exp(xi*(t-t1-nu))));
  }
  
  real[] SEIR(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    
    int K = x_i[1];
    real tmax2 = x_r[1];
    real tswitch = x_r[2];
    real dydt[(6*K)];                //SEIAR, then C 
    real nI = sum(y[(2*K+1):(3*K)]); // total infectious
    real ntot = sum(y[1:(4*K)]);
    
    real beta = theta[1];  // transmission rate
    real eta = theta[2];   // reduction in transmission rate after quarantine
    real xi = theta[3];    // slope of quarantine implementation
    real nu = theta[4];    // shift of quarantine implementation
    real tau = 1.0/x_r[3]; // incubation period
    real mu = 1.0/x_r[4];  // infectious period
    real psi = x_r[5];     // probability of symptoms

    real p_tmax = switch_zero(t,tmax2); // tau is 0 after tmax (stop "recruiting")
    real p_tswitch = switch_eta(t,tswitch,eta,nu,xi);
    
    for (k in 1:K) {
      dydt[k] = -beta * p_tswitch * y[k] * nI/ntot; // S
      dydt[K+k] = beta * p_tswitch * y[k] * nI/ntot - tau * p_tmax * y[K+k]; // E
      dydt[2*K+k] = psi * tau * p_tmax * y[K+k] - mu * y[2*K+k]; // I
      dydt[3*K+k] = (1-psi) * tau * p_tmax * y[K+k] - mu * y[3*K+k]; // A
      dydt[4*K+k] =  mu * y[2*K+k] +  mu * y[3*K+k]; // R
      dydt[5*K+k] = psi * tau * p_tmax * y[K+k];  // C
    }
    return(dydt);
    
  }

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

  // Replacement for ode_integrate_rk45
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
  int K;               //number of age classes
  vector[K] age_dist;
  int pop_t;           //total population
  int tmax;            //number of days 
  real tmax2;
  real tswitch;
  
  // Data to fit
  int D;                   //number of days with reported incidence
  int incidence_cases[D];  //overal incidence for W weeks
  int incidence_deaths[D]; //overal incidence for W weeks
  int agedistr_cases[K];   //number of cases at tmax for the K age classes
  int agedistr_deaths[K];  //mortality at tmax for the K age classes
  
  // Parameters in priors
  real p_beta;
  real p_eta[2];
  real p_pi[2];
  real p_epsilon[2];
  real p_rho[2];
  real p_phi;
  real p_xi;
  real p_nu;
  
  // Fixed parameters
  real p_incubation;
  real p_infectious;
  real p_psi;
  int G;
  real p_gamma[G];
  
  // Simulation
  real t0;    //starting time
  int t_data; //time of first data
  int S;
  real ts[S]; // time bins
  int inference;
  int doprint;
}

transformed data {
  real x_r[5] = {tmax2, tswitch, p_incubation, p_infectious, p_psi};
  int x_i[1] = {K};
  real EPS = 1.0E-9;
}

parameters{
  real<lower=0> beta;                  // base transmission rate
  real<lower=0,upper=1> eta;           // reduct. in transm. rate after quarant. measures
  vector<lower=0,upper=1> [K] epsilon; // age-dependent mortality probability
  vector<lower=0,upper=1> [K-1] rho;   // age-dependent reporting probability
  real<lower=0, upper=1> pi;           // number of cases at t0
  real<lower=0> phi[2];                // variance parameters
  real<lower=0,upper=1> xi_raw;        // slope of quarantine implementation
  real<lower=0> nu;                    // shift of quarantine implementation
}

transformed parameters {
  vector[K] rho_K;
  real theta[5];         // ODE parameters
  real xi = xi_raw+0.5;
  for(i in 1:(K-1)){ 
    rho_K[i] = rho[i];
  }
  rho_K[K] = 1.0;
  theta[1:5] = {beta, eta, xi, nu, 0.0};
}

model {
  real log_lh;
#include inference_midpoint.stan
  
  // priors
  target += exponential_lpdf(beta | p_beta);
  target += exponential_lpdf(nu | p_nu);
  target += exponential_lpdf(phi | p_phi);
  target += beta_lpdf(eta | p_eta[1], p_eta[2]);
  target += beta_lpdf(pi | p_pi[1], p_pi[2]);
  target += beta_lpdf(xi_raw | p_xi, p_xi);
  for(k in 1:K) { target += beta_lpdf(epsilon[k] | p_epsilon[1], p_epsilon[2]); }
  for(k in 1:(K-1)){ target += beta_lpdf(rho[k] | p_rho[1], p_rho[2]); }
  
  // likelihood
#include log_lh.stan
  target += log_lh;
}

generated quantities{
  real logp_1;
  real logp_2;
  {
    real log_lh = 0.0;
#include inference_bdf.stan
    logp_1 = log_lh;
  }
  {
    real log_lh = 0.0;
#include inference_midpoint.stan
    logp_2 = log_lh;
  }
}
