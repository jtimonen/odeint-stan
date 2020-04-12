functions {
#include stan_functions/switch.stan
#include stan_functions/SEIR.stan
#include stan_functions/extract.stan
#include stan_functions/likelihood.stan
#include stan_functions/prior.stan
}

data {
#include stan_chunks/data.stan
}

transformed data {
#include stan_chunks/transformed_data.stan
}

parameters{
#include stan_chunks/parameters.stan
}

transformed parameters {
#include stan_chunks/transformed_parameters.stan
  y = integrate_ode_bdf(SEIR, init, t0, ts, theta, x_r, x_i, abs_tol, rel_tol, max_iter);
#include stan_chunks/prior_and_likelihood.stan
}

model {
  // evaluate lp__ with the (invisible) jacobian adjustment term included
  target += log_prior_na;
  target += log_lik_na;
}

generated quantities{
  
}
