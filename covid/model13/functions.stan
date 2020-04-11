functions {
  
  // Sigmoidal switch
  real switch_eta(real t, real t1, real eta, real nu, real xi) {
    return(eta+(1-eta)/(1+exp(xi*(t-t1-nu))));
  }
  
  // Prior
  real log_prior_noadjustment(
      real beta, 
      real eta, 
      vector epsilon, 
      vector raw_rho,
      real xi_raw, 
      real pii, 
      real psi, 
      real[] phi, 
      real nu,
      data real p_beta,
      data real[] p_eta,
      data real[] p_epsilon,
      data real[] p_rho,
      data real p_xi,
      data real[] p_pi,
      data real[] p_psi,
      data real p_phi,
      data real p_nu)
  {
    real log_prior = 0.0;
    int K = num_elements(epsilon);
    log_prior += beta_lpdf(beta | p_beta, p_beta);
    log_prior += beta_lpdf(eta | p_eta[1], p_eta[2]);
    for(k in 1:K){ 
      log_prior += beta_lpdf(epsilon[k] | p_epsilon[1], p_epsilon[2]);
    }
    for(k in 1:(K-1)) {
      log_prior += beta_lpdf(raw_rho[k] | p_rho[1], p_rho[2]);
    }
    log_prior += beta_lpdf(xi_raw | p_xi, p_xi); 
    log_prior += beta_lpdf(pii | p_pi[1], p_pi[2]);
    log_prior += beta_lpdf(psi | p_psi[1], p_psi[2]);
    for(j in 1:2){
      log_prior += exponential_lpdf(phi[j] | p_phi);
    }
    log_prior += exponential_lpdf(nu | p_nu);
    return(log_prior);
  }
  
  // Likelihood
  real log_likelihood_noadjustment(
    real psi, 
    real pii, 
    real[] phi, 
    real[] p_gamma, 
    vector epsilon,  
    vector rho,  
    data int[] icases, 
    data int[] ideaths, 
    data int[] acases, 
    data int[] adeaths,
    data vector age_dist, 
    data real EPS, 
    data int pop_t, 
    data int t_data, 
    real[,] y){
    
    // Declarations
    real log_lik = 0.0;
    
    int S = size(y);
    int G = num_elements(p_gamma);
    int K = num_elements(age_dist);
    int D = num_elements(icases);
    
    vector[K] comp_S[S];
    vector[K] comp_E[S];
    vector[K] comp_I[S];
    vector[K] comp_C[S+G];
    vector[K] comp_diffC[S+G];
    vector[K] comp_diffM[S+G];
    vector[K] comp_M[S+G];
    
    vector[D] out_icases;   // overall case incidence by day
    vector[D] out_ideaths;  // overal mortality incidence by day 
    vector[K] out_acases;   // final age distribution of cases (simplex)
    vector[K] out_adeaths;  // final age distribution of deaths (simplex)
    
    // extract and format ODE results 
    for(i in 1:S) {
      comp_S[i] = (to_vector(y[i,1:K]) + to_vector(age_dist) * (1-pii) + EPS) * pop_t;
      comp_E[i] = (to_vector(y[i,(K+1):(2*K)]) + to_vector(age_dist) * pii + EPS) * pop_t;
      comp_I[i] = (to_vector(y[i,(2*K+1):(3*K)]) + EPS) * pop_t;
      comp_C[i] = (to_vector(y[i,(3*K+1):(4*K)]) + EPS) * pop_t;
      comp_diffC[i] = i==1 ? comp_C[i,] : EPS*pop_t + comp_C[i,] - comp_C[i-1,]; 
      // comp_diffC = lagged difference of cumulative incidence of symptomatics
    }
    
    // Incidence and cumulative incidence after S
    for(g in 1:G){
      comp_C[S+g] = comp_C[S];
      comp_diffC[S+g] = rep_vector(EPS, K);
    }
    
    // Mortality, set diffM and M to 0
    for(i in 1:(G+S)){
      comp_diffM[i] = rep_vector(EPS, K);
    }
    
    // Compute mortality
    for(i in 1:S) {
      for(g in 1:G) {
        comp_diffM[i+g] += comp_diffC[i] .* epsilon * p_gamma[g] ;
      }
    }
    
    // Cumulative sum
    for(i in 1:(S+G)) {
      for(k in 1:K) {
        comp_M[i,k] = sum(comp_diffM[1:i,k]);
      }
    }
    
    // Compute outcomes
    for(i in t_data:S){
      out_icases[i-t_data+1] = sum(comp_diffC[i].*rho);
      out_ideaths[i-t_data+1] = sum(comp_diffM[i]);
    }
    out_acases = (comp_C[S,].*rho) ./ sum(comp_C[S,].*rho);
    out_adeaths = (comp_M[S,]) ./ sum(comp_M[S,]);
    
    // Evaluate likelihood
    for(i in 1:D) {
      log_lik += neg_binomial_2_lpmf( icases[i] | out_icases[i], out_icases[i]/phi[1]);
      log_lik  += neg_binomial_2_lpmf( ideaths[i] | out_ideaths[i], out_ideaths[i]/phi[2]);
    }
    log_lik += multinomial_lpmf(acases | out_acases);
    log_lik += multinomial_lpmf(adeaths | out_adeaths);
    return log_lik;
  }
  
  // ODE system
  real[] SEIR(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    int K = x_i[1];
    real tswitch = x_r[1];
    real dydt[(4*K)];  // SEIAR, then C and D
    real nI;           // total infectious
    real ntot;
    real beta;         // transmission rate
    real eta;          // reduction in transmission rate after quarantine
    real xi;           // slope of quarantine implementation
    real nu;           // shift of quarantine implementation
    real tau;          // incubation period
    real mu;           // infectious period
    real psi;          // probability of symptoms
    real p_tswitch;
    real contact[K*K]; // contact matrix, first K values
    // corresponds to number of contacts between 
    // age class 1 and other classes, etc
    
    real n_by_age[K];
    real f_inf[K];     // force of infection
    real init[K*4];
    real age_dist[K];
    real pii;           // number of cases at t0
    
    // estimated parameters
    beta = theta[1];
    eta = theta[2];
    xi = theta[3];
    nu = theta[4];
    pii = theta[5];
    psi = theta[6];
    
    // Initial conditions
    for(k in 1:K){
      age_dist[k] = x_r[3+K*K + k];
      init[k] = age_dist[k] * (1-pii);
      init[K+k] = age_dist[k] * pii;
      init[2*K+k] = 0.0;
      init[3*K+k] = 0.0;
    }
    
    // Fixed parameters
    tau = 1.0/x_r[2];
    mu = 1.0/x_r[3];
    contact = x_r[4:(3+K*K)];
    
    // Total number of infectious people
    p_tswitch = switch_eta(t, tswitch, eta, nu, xi);
    
    // Force of infection by age classes: 
    // beta * p_tswitch * sum((#infected) / (#people) * (#contacts))
        for(k in 1:K){
          f_inf[k] = beta * p_tswitch * sum(to_vector(y[(2*K+1):(3*K)]) ./ to_vector(age_dist) .* to_vector(contact[(K*(k-1)+1):(k*K)])); //
        }
        for (k in 1:K) {
          dydt[k] = - f_inf[k] * (y[k]+init[k]);                                      // S
          dydt[K+k] = f_inf[k] * (y[k]+init[k])- tau * (y[K+k]+init[K+k]);            // E
          dydt[2*K+k] = psi * tau * (y[K+k]+init[K+k]) - mu * (y[2*K+k]+init[2*K+k]); // I
          dydt[3*K+k] = psi * tau * (y[K+k]+init[K+k]);                               // C
        }
        return(dydt);
  }
}
