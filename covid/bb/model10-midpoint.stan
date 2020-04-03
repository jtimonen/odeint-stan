functions {
  real switch_zero(real t, real t1) {
    return(1/(1+exp(5*(t-t1))));
  }
  real switch_eta(real t, real t1, real eta, real nu, real xi) {
    return(eta+(1-eta)/(1+exp(xi*(t-t1-nu))));
  }
  real[] SEIR(real t,
  real[] y,
  real[] theta,
  real[] x_r,
  int[] x_i
  ) {
    int K = x_i[1];
    real tmax2 = x_r[1];
    real tswitch = x_r[2];
    real dydt[(6*K)]; //SEIAR, then C 
    real nI; // total infectious
    real ntot;
    
    real beta; // transmission rate
    real eta; // reduction in transmission rate after quarantine
    real xi; // slope of quarantine implementation
    real nu; // shift of quarantine implementation
    real tau; // incubation period
    real mu; // infectious period
    real psi; // probability of symptoms
    real p_tmax; // tau is 0 after tmax (stop "recruiting")
    real p_tswitch;
    
    // estimated parameters
    beta = theta[1];
    eta = theta[2];
    xi = theta[3];
    nu = theta[4];
    
    // fixed parameters
    tau = 1.0/x_r[3];
    mu = 1.0/x_r[4];
    psi = x_r[5];
    
    //Total number of infectious people
    nI = sum(y[(2*K+1):(3*K)]);
    ntot = sum(y[1:(4*K)]);
    p_tmax = switch_zero(t,tmax2);
    p_tswitch = switch_eta(t,tswitch,eta,nu,xi);
    
    for (k in 1:K) {
      // S
      dydt[k] = -beta * p_tswitch * y[k] * nI/ntot; 
      // E
      dydt[K+k] = beta * p_tswitch * y[k] * nI/ntot - tau * p_tmax * y[K+k];
      // I
      dydt[2*K+k] = psi * tau * p_tmax * y[K+k] - mu * y[2*K+k];
      // A
      dydt[3*K+k] = (1-psi) * tau * p_tmax * y[K+k] - mu * y[3*K+k];
      // R
      dydt[4*K+k] =  mu * y[2*K+k] +  mu * y[3*K+k];
      // C
      dydt[5*K+k] = psi * tau * p_tmax * y[K+k];
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
  int K; //number of age classes
  vector[K] age_dist;
  int pop_t; //total population
  int tmax; //number of days between date of first "infection" (late November) and 11 February (when data is collected)
  real tmax2;
  real tswitch;
  
  //Data to fit
  int D; //number of days with reported incidence
  int incidence_cases[D]; //overal incidence for W weeks
  int incidence_deaths[D]; //overal incidence for W weeks
  int agedistr_cases[K]; //number of cases at tmax for the K age classes
  int agedistr_deaths[K]; //mortality at tmax for the K age classes
  //Parameters in priors
  real p_beta;
  real p_eta[2];
  real p_pi[2];
  real p_epsilon[2];
  real p_rho[2];
  real p_phi;
  real p_xi;
  real p_nu;
  // fixed parameters
  real p_incubation;
  real p_infectious;
  real p_psi;
  int G;
  real p_gamma[G];
  //Simulation
  real t0; //starting time
  int t_data; //time of first data
  int S;
  real ts[S]; // time bins
  
  int inference;
  int doprint;
}

transformed data {
  real x_r[5] = {tmax2,tswitch,p_incubation,p_infectious,p_psi};
  int x_i[1] = {K};
}

parameters{
  real<lower=0> beta; // base transmission rate
  real<lower=0,upper=1> eta; // reduction in transmission rate after quarantine measures
  vector<lower=0,upper=1> [K] epsilon; // age-dependent mortality probability
  vector<lower=0,upper=1> [K-1] rho; // age-dependent reporting probability
  real<lower=0, upper=1> pi; // number of cases at t0
  real<lower=0> phi[2]; // variance parameters
  real<lower=0,upper=1> xi_raw; // slope of quarantine implementation
  real<lower=0> nu; // shift of quarantine implementation
}
transformed parameters {
  // transformed parameters
  vector[K] rho_K;
  real xi = xi_raw+0.5;
  // change of format for integrate_ode_rk45
  real theta[5]; // vector of parameters
  real y[S,K*6]; // raw ODE output
  vector[K] comp_S[S];
  vector[K] comp_E[S];
  vector[K] comp_I[S];
  vector[K] comp_A[S];
  vector[K] comp_R[S];
  vector[K] comp_C[S];
  vector[K] comp_diffC[S];
  vector[K] comp_diffM[S];
  vector[K] comp_M[S];
  // outcomes
  vector[D] output_incidence_cases; // overall case incidence by day
  vector[D] output_incidence_deaths; // overal mortality incidence by day 
  simplex[K] output_agedistr_cases; // final age distribution of cases
  simplex[K] output_agedistr_deaths; // final age distribution of deaths
  
  // transformed parameters
  for(i in 1:(K-1)){
    rho_K[i]=rho[i];
  }
  rho_K[K]=1.0;
  // change of format for integrate_ode_rk45
  theta[1:5] = {beta,eta,xi,nu, 0.0};
    {
  real init[K*6]; // initial values
  for(k in 1:K){
    init[k] = age_dist[k] * (1-pi);
    init[K+k] = age_dist[k] * pi;
    init[2*K+k] = 0.0;
    init[3*K+k] = 0.0;
    init[4*K+k] = 0.0;
    init[5*K+k] = 0.0;
  }
  // run ODE solver
  /*y = integrate_ode_bdf(
    SEIR, // ODE function
    init, // initial states
    t0, // t0
    ts, // evaluation dates (ts)
    theta, // parameters
    x_r, // real data
    x_i, // integer data
    1.0E-6, 1.0E-6, 1.0E3);*/ // tolerances and maximum steps
  // Comment out the previous integrator and use this
  //  to use the explicit midpoint integrator.
  //  The third to last argument is the timestep
  y = explicit_midpoint_integrator(
    init, // initial states
    t0, // t0
    ts, // evaluation dates (ts)
    theta, // parameters
    5e-1, // stepsize
    x_r, // real data
    x_i); // tolerances and maximum steps
    }
    // extract and format ODE results (1.0E-9 correction to avoid negative values due to unprecise estimates of zeros as tolerance is 1.0E-10)
    for(i in 1:S) {
      comp_S[i] = (to_vector(y[i,1:K]) + 1.0E-9) * pop_t;
      comp_E[i] = (to_vector(y[i,(K+1):(2*K)]) + 1.0E-9) * pop_t;
      comp_I[i] = (to_vector(y[i,(2*K+1):(3*K)]) + 1.0E-9) * pop_t;
      comp_A[i] = (to_vector(y[i,(3*K+1):(4*K)]) + 1.0E-9) * pop_t;
      comp_R[i] = (to_vector(y[i,(4*K+1):(5*K)]) + 1.0E-9) * pop_t;
      comp_C[i] = (to_vector(y[i,(5*K+1):(6*K)]) + 1.0E-9) * pop_t;
      comp_diffC[i] = i==1 ? comp_C[i,] : 1.0E-9*pop_t + comp_C[i,] - comp_C[i-1,]; // lagged difference of cumulative incidence
    }
    // compute mortality
    for(i in 1:S) {
      comp_diffM[i] = rep_vector(1.0E-9,K);
      comp_M[i] = rep_vector(0.0,K);
      for(g in 1:G) {
        if (i>g) comp_diffM[i] += comp_diffC[i-g] .* epsilon * p_gamma[g] ;
      }
      comp_M[i] += comp_diffM[i];
    }
    // compute outcomes (again, 1.0E-9 correction to avoid negative values in the lagged differences)
    for(i in t_data:tmax){
      output_incidence_cases[i-t_data+1] = sum(comp_diffC[i].*rho_K);
      output_incidence_deaths[i-t_data+1] = sum(comp_diffM[i]);
    }
    output_agedistr_cases = (comp_C[tmax,].*rho_K) ./ sum(comp_C[tmax,].*rho_K);
    output_agedistr_deaths = (comp_M[tmax,]) ./ sum(comp_M[tmax,]);
}
model {
  // priors
  beta ~ exponential(p_beta);
  eta ~ beta(p_eta[1],p_eta[2]);
  for(k in 1:K) epsilon[k] ~ beta(p_epsilon[1],p_epsilon[2]);
  for(k in 1:(K-1)) rho[k] ~ beta(p_rho[1],p_rho[2]);
  pi ~ beta(p_pi[1],p_pi[2]);
  phi ~ exponential(p_phi);
  xi_raw ~ beta(p_xi,p_xi); 
  nu ~ exponential(p_nu);
  // debug
  if(doprint==1) {
    print("beta: ",beta);
    print("eta: ",beta);
    print("epsilon: ",epsilon);
    print("rho: ",rho);
    print("pi: ",pi);
    print("y[5,]: ",y[5,]);
    print("comp_C[5,]: ",comp_C[5,]);
    print("comp_diffC[5,]: ",comp_diffC[5,]);
    print("comp_M[5,]: ",comp_M[5,]);
    print("comp_diffM[5,]: ",comp_diffM[5,]);
  }
  // likelihood
  if (inference==1) {
    for(i in 1:D) {
      target += neg_binomial_2_lpmf( incidence_cases[i] | output_incidence_cases[i], output_incidence_cases[i]/phi[1]);
      target += neg_binomial_2_lpmf( incidence_deaths[i] | output_incidence_deaths[i],output_incidence_deaths[i]/phi[2]);
    }
    target += multinomial_lpmf(agedistr_cases | output_agedistr_cases);
    target += multinomial_lpmf(agedistr_deaths | output_agedistr_deaths);
  }
}

/*generated quantities{
  real max_abs_error = 0.0;
  real max_rel_error = 0.0;
  
  {
    real y_ref[S,K*6];
  
    y_ref = integrate_ode_bdf(
			      SEIR, // ODE function
			      init, // initial states
			      t0, // t0
			      ts, // evaluation dates (ts)
			      theta, // parameters
			      x_r, // real data
			      x_i, // integer data
			      1.0E-10, 1.0E-10, 1.0E3);

    for(s in 1:S) {
      for(k in 1:(6 * K)) {
	max_abs_error = max([max_abs_error, abs(y[s, k] - y_ref[s, k])]);
	max_rel_error = max([max_rel_error, abs(y[s, k] - y_ref[s, k]) /
			     max([abs(y[s, k]), abs(y_ref[s, k])])]);
      }
    }
  }
}
*/
