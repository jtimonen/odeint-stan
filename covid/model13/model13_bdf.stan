functions {
  real switch_eta(real t, real t1, real eta, real nu, real xi) {
    return(eta+(1-eta)/(1+exp(xi*(t-t1-nu))));
  }
  real[] SEIR(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    int K = x_i[1];
    real tswitch = x_r[1];
    real dydt[(4*K)]; //SEIAR, then C and D
    real nI; // total infectious
    real ntot;
    
    real beta; // transmission rate
    real eta;  // reduction in transmission rate after quarantine
    real xi;   // slope of quarantine implementation
    real nu;   // shift of quarantine implementation
    real tau;  // incubation period
    real mu;   // infectious period
    real psi;  // probability of symptoms
    real p_tswitch;
    real contact[K*K]; //contact matrix, first K values, corresponds to number of contact between age class 1 and other classes, etc
    real n_by_age[K];
    real f_inf[K];
    
    real init[K*4];
    real age_dist[K];
    real pi; // number of cases at t0
    
    // estimated parameters
    beta = theta[1];
    eta = theta[2];
    xi = theta[3];
    nu = theta[4];
    pi = theta[5];
    psi = theta[6];
    
    // Initial conditions
    for(k in 1:K){
      age_dist[k] = x_r[3+K*K + k];
      init[k] = age_dist[k] * (1-pi);
      init[K+k] = age_dist[k] * pi;
      init[2*K+k] = 0.0;
      init[3*K+k] = 0.0;
    }
    
    // fixed parameters
    tau = 1.0/x_r[2];
    mu = 1.0/x_r[3];
    contact = x_r[4:(3+K*K)];
    
    //Total number of infectious people
    p_tswitch = switch_eta(t,tswitch,eta,nu,xi);
    
    //Force of infection by age classes: beta * p_tswitch * sum((number of infected people by age) / (total number of people by age) * (number of contact by age))
    for(k in 1:K){
      f_inf[k] = beta * p_tswitch * sum(to_vector(y[(2*K+1):(3*K)]) ./ to_vector(age_dist) .* to_vector(contact[(K*(k-1)+1):(k*K)])); //
    }
    for (k in 1:K) {
      // S
      dydt[k] = - f_inf[k] * (y[k]+init[k]); 
      // E
      dydt[K+k] = f_inf[k] * (y[k]+init[k])- tau * (y[K+k]+init[K+k]);
      // I
      dydt[2*K+k] = psi * tau * (y[K+k]+init[K+k]) - mu * (y[2*K+k]+init[2*K+k]);
      // C
      dydt[3*K+k] = psi * tau * (y[K+k]+init[K+k]);
    }
    return(dydt);
  }
}

data {
  int K; //number of age classes
  vector[K] age_dist;
  int pop_t; //total population
  real tswitch;
  
  //Data to fit
  int D;                   // number of days with reported incidence
  int incidence_cases[D];  // overal incidence for W weeks
  int incidence_deaths[D]; // overal incidence for W weeks
  int agedistr_cases[K];   // number of cases at tmax for the K age classes
  int agedistr_deaths[K];  // mortality at tmax for the K age classes
  //Parameters in priors
  real p_beta;
  real p_eta[2];
  real p_pi[2];
  real p_epsilon[2];
  real p_rho[2];
  real p_phi;
  real p_xi;
  real p_nu;
  real p_psi[2];
  // fixed parameters
  real p_incubation;
  real p_infectious;
  int G;
  real p_gamma[G];
  //Simulation
  real t0; //starting time
  int t_data; //time of first data
  int S;
  real ts[S]; // time bins
  real contact[K*K];
  
  int inference;
  int doprint;
}

transformed data {
  real EPS = 1.0E-9; // correction to avoid neg. values as tolerance is 1.0E-10
  real x_r[3+K*K+K]; //4 parameters + K*K contact matrix parameters + K age_dist parameters
  int x_i[1] = {K};
  real init[K*4] = rep_array(0.0, K * 4); // initial values
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
  real<lower=0, upper=1> pi;              // number of cases at t0
  real<lower=0> phi[2];                   // variance parameters
  real<lower=0,upper=1> xi_raw;           // slope of quarantine implementation
  real<lower=0> nu;                       // shift of quarantine implementation
  real<lower=0,upper=1> psi;              // proportion of symptomatics
}

transformed parameters {
  // transformed parameters
  vector[K] rho;
  real xi = xi_raw+0.5;
  // change of format for integrate_ode_rk45
  real theta[6]; // vector of parameters
  real y[S,K*4]; // raw ODE output
  vector[K] comp_S[S];
  vector[K] comp_E[S];
  vector[K] comp_I[S];
  vector[K] comp_D[S];
  vector[K] comp_diffD[S];
  vector[K] comp_C[S+G];
  vector[K] comp_diffC[S+G];
  vector[K] comp_diffM[S+G];
  vector[K] comp_M[S+G];
  // outcomes
  vector[D] output_incidence_cases;  // overall case incidence by day
  vector[D] output_incidence_deaths; // overal mortality incidence by day 
  simplex[K] output_agedistr_cases;  // final age distribution of cases
  simplex[K] output_agedistr_deaths; // final age distribution of deaths
  
  // transformed paremeters
  for(i in 1:(K-1)){
    rho[i]=raw_rho[i];
  }
  rho[K] = 1.0;
  // change of format for integrate_ode_rk45
  theta[1:6] = {beta,eta,xi,nu,pi,psi};
  // run ODE solver
  y = integrate_ode_bdf(SEIR, init, t0, ts, theta, x_r, x_i, 1.0E-10, 1.0E-10, 1.0E3);
    // extract and format ODE results 
    for(i in 1:S) {
      comp_S[i] = (to_vector(y[i,1:K]) + to_vector(age_dist) * (1-pi) + EPS) * pop_t;
      comp_E[i] = (to_vector(y[i,(K+1):(2*K)]) + to_vector(age_dist) * pi + EPS) * pop_t;
      comp_I[i] = (to_vector(y[i,(2*K+1):(3*K)]) + EPS) * pop_t;
      comp_C[i] = (to_vector(y[i,(3*K+1):(4*K)]) + EPS) * pop_t;
      comp_diffC[i] = i==1 ? comp_C[i,] : EPS*pop_t + comp_C[i,] - comp_C[i-1,]; // lagged difference of cumulative incidence of symptomatics
    }
    //Incidence and cumulative incidence after S
    for(g in 1:G){
      comp_C[S+g] = comp_C[S];
      comp_diffC[S+g] = rep_vector(EPS, K);
    }
    //Mortality
    //set diffM and M to 0
    for(i in 1:(G+S)){
      comp_diffM[i] = rep_vector(EPS, K);
    }
    //compute mortality
    for(i in 1:S) {
      for(g in 1:G) {
        comp_diffM[i+g] += comp_diffC[i] .* epsilon * p_gamma[g] ;
      }
    }
    // cumulative sum
    for(i in 1:(S+G)) {
      for(k in 1:K) {
        comp_M[i,k] = sum(comp_diffM[1:i,k]);
      }
    }
    //compute D and diffD
    for(i in 1:S){
      comp_D[i] = (1.0-psi)/psi * comp_C[i];
      comp_diffD[i] = (1.0-psi)/psi * comp_diffC[i];
    }
    // compute outcomes
    for(i in t_data:S){
      output_incidence_cases[i-t_data+1] = sum(comp_diffC[i].*rho);
      output_incidence_deaths[i-t_data+1] = sum(comp_diffM[i]);
    }
    output_agedistr_cases = (comp_C[S,].*rho) ./ sum(comp_C[S,].*rho);
    output_agedistr_deaths = (comp_M[S,]) ./ sum(comp_M[S,]);
}

model {
  // priors
  beta ~ beta(p_beta,p_beta);
  eta ~ beta(p_eta[1],p_eta[2]);
  for(k in 1:K) epsilon[k] ~ beta(p_epsilon[1],p_epsilon[2]);
  for(k in 1:(K-1)) raw_rho[k] ~ beta(p_rho[1],p_rho[2]);
  pi ~ beta(p_pi[1],p_pi[2]);
  phi ~ exponential(p_phi);
  xi_raw ~ beta(p_xi,p_xi); 
  nu ~ exponential(p_nu);
  psi ~ beta(p_psi[1],p_psi[2]);

  // likelihood
  for(i in 1:D) {
      target += neg_binomial_2_lpmf( incidence_cases[i] | output_incidence_cases[i], output_incidence_cases[i]/phi[1]);
      target += neg_binomial_2_lpmf( incidence_deaths[i] | output_incidence_deaths[i],output_incidence_deaths[i]/phi[2]);
  }
  target += multinomial_lpmf(agedistr_cases | output_agedistr_cases);
  target += multinomial_lpmf(agedistr_deaths | output_agedistr_deaths);
}
