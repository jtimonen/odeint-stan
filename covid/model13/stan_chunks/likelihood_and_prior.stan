  // Compute prior and likelihood
  log_prior_na += log_prior_noadjustment(beta, eta, epsilon, raw_rho, xi_raw, 
      pii, psi, phi, nu, p_beta, p_eta, p_epsilon, p_rho, p_xi, p_pi, 
      p_psi, p_phi, p_nu);
  log_lik_na += log_likelihood_noadjustment(psi, pii, phi, p_gamma, epsilon, raw_rho,
      incidence_cases, incidence_deaths, agedistr_cases, agedistr_deaths, age_dist,
      EPS, pop_t, t_data, y);
