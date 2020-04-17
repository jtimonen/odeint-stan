  real log_prior_na_REF_ = 0.0;
  real log_lik_na_REF_ = 0.0;
  
  // Solve ODE
  real y_REF_[S,K*4] = integrate_ode_bdf(SEIR, init, t0, ts, theta, x_r, x_i,
      abs_tol_REF_, rel_tol_REF_, max_iter_REF_);
  
  // Compute prior and likelihood
  log_prior_na_REF_ += log_prior_noadjustment(beta, eta, epsilon, raw_rho, 
      xi_raw, pii, psi, phi, nu, p_beta, p_eta, p_epsilon, p_rho, p_xi, p_pi,
      p_psi, p_phi, p_nu);
  log_lik_na_REF_ += log_likelihood_noadjustment(psi, pii, phi, p_gamma, epsilon, raw_rho,
      incidence_cases, incidence_deaths, agedistr_cases, agedistr_deaths, age_dist,
      EPS, pop_t, t_data, y_REF_);
