  real y[S,K*6];         // raw ODE output
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
  vector[D] output_incidence_cases;  // overall case incidence by day
  vector[D] output_incidence_deaths; // overal mortality incidence by day 
  vector[K] output_agedistr_cases;  // final age distribution of cases
  vector[K] output_agedistr_deaths; // final age distribution of deaths
  
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
    y = integrate_ode_bdf(SEIR, init, t0, ts, theta, x_r, x_i, 1.0E-6, 1.0E-6, 1.0E3);
    // y = explicit_midpoint_integrator(init, t0, ts, theta, 5e-1, x_r, x_i);
  }
  
  // extract and format ODE results
  for(i in 1:S) {
    comp_S[i] = (to_vector(y[i,1:K]) + EPS) * pop_t;
    comp_E[i] = (to_vector(y[i,(K+1):(2*K)]) + EPS) * pop_t;
    comp_I[i] = (to_vector(y[i,(2*K+1):(3*K)]) + EPS) * pop_t;
    comp_A[i] = (to_vector(y[i,(3*K+1):(4*K)]) + EPS) * pop_t;
    comp_R[i] = (to_vector(y[i,(4*K+1):(5*K)]) + EPS) * pop_t;
    comp_C[i] = (to_vector(y[i,(5*K+1):(6*K)]) + EPS) * pop_t;
    comp_diffC[i] = i==1 ? comp_C[i,] : EPS*pop_t + comp_C[i,] - comp_C[i-1,];
  }
  
  // compute mortality
  for(i in 1:S) {
    comp_diffM[i] = rep_vector(EPS, K);
    comp_M[i] = rep_vector(0.0,K);
    for(g in 1:G) {
      if (i>g) comp_diffM[i] += comp_diffC[i-g] .* epsilon * p_gamma[g] ;
    }
    comp_M[i] += comp_diffM[i];
  }
  
  // compute outcomes
  for(i in t_data:tmax){
    output_incidence_cases[i-t_data+1] = sum(comp_diffC[i].*rho_K);
    output_incidence_deaths[i-t_data+1] = sum(comp_diffM[i]);
  }
  output_agedistr_cases = (comp_C[tmax,].*rho_K) ./ sum(comp_C[tmax,].*rho_K);
  output_agedistr_deaths = (comp_M[tmax,]) ./ sum(comp_M[tmax,]);
