
// Likelihood
real log_likelihood_noadjustment(
    real psi, 
    real pii, 
    real[] phi, 
    real[] p_gamma, 
    vector epsilon,  
    vector raw_rho,  
    data int[] icases, 
    data int[] ideaths, 
    data int[] acases, 
    data int[] adeaths,
    data vector age_dist, 
    data real EPS, 
    data int pop_t, 
    data int t_data, 
    real[,] y)
{
  
  // Declare variables
  real log_lik = 0.0;
  int S = size(y);
  int G = num_elements(p_gamma);
  int K = num_elements(age_dist);
  int D = num_elements(icases);
  vector[D] out_icases;   // overall case incidence by day
  vector[D] out_ideaths;  // overal mortality incidence by day 
  vector[K] out_acases;   // final age distribution of cases (simplex)
  vector[K] out_adeaths;  // final age distribution of deaths (simplex)
  
  // Extract compartments from y
  vector[K] comp_S[S] = get_comp_S(y, age_dist, pii, pop_t, EPS);
  vector[K] comp_E[S] = get_comp_E(y, age_dist, pii, pop_t, EPS);
  vector[K] comp_I[S] = get_comp_I(y, pop_t, EPS, K);
  vector[K] comp_C[S+G] = get_comp_C(y, pop_t, EPS, K, G);
  vector[K] comp_diffC[S+G] = get_comp_diffC(comp_C, pop_t, EPS, S, G);
  vector[K] comp_diffM[S+G] = get_comp_diffM(comp_diffC, epsilon, p_gamma, EPS, S, G);
  vector[K] comp_M[S+G] = get_comp_M(comp_diffM, S, G);

  // Append rho of last age class
  vector[K] rho;
  for(i in 1:(K-1)){
    rho[i] = raw_rho[i];
  }
  rho[K] = 1.0;
  
  // Compute outcomes 
  for(i in t_data:S){
    out_icases[i-t_data+1] = sum(comp_diffC[i].*rho);
    out_ideaths[i-t_data+1] = sum(comp_diffM[i]);
  }
  out_acases = (comp_C[S,].*rho) ./ sum(comp_C[S,].*rho);
  out_adeaths = (comp_M[S,]) ./ sum(comp_M[S,]);

  // Evaluate likelihood
  for(i in 1:D) {
    log_lik += neg_binomial_2_lpmf( icases[i] | out_icases[i], out_icases[i] / phi[1] );
    log_lik  += neg_binomial_2_lpmf( ideaths[i] | out_ideaths[i], out_ideaths[i] / phi[2] );
  }
  log_lik += multinomial_lpmf(acases | out_acases);
  log_lik += multinomial_lpmf(adeaths | out_adeaths);
  return(log_lik);
}
