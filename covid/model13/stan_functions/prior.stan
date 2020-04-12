
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
