functions {
#include stan_functions/switch.stan
#include stan_functions/ODE.stan
#include stan_functions/extract.stan
#include stan_functions/likelihood.stan
#include stan_functions/prior.stan

  // Explicit midpoint method
  real[,] integrate_ode_emp(real[] y0, real t0, real[] ts, real[] theta,
      data real STEP_SIZE, data int[] INTERP_R, data real[] INTERP_A, 
      data real[] x_r, data int[] x_i){
    int d = size(y0);
    int n = size(ts);
    real x[n, d];
    int R_n = INTERP_R[n];
    real t_tmp = t0;
    real t_mid;
    vector[d] y_tmp = to_vector(y0);
    vector[d] y_mid;
    vector[d] y[R_n+1];
    y[1] = y_tmp;
    for(i in 1:R_n){
      y_mid = y_tmp + 0.5*STEP_SIZE*odefun(t_tmp, y_tmp, theta, x_r, x_i);
      t_mid = t_tmp + 0.5*STEP_SIZE;
      y_tmp = y_tmp + STEP_SIZE*odefun(t_mid, y_mid, theta, x_r, x_i);
      t_tmp = t_tmp + STEP_SIZE;
      y[i+1] = y_tmp;
    }
    x = interpolate(y, INTERP_R, INTERP_A);
    return(x);
  }
  
}

data {
#include stan_chunks/data.stan
#include stan_chunks/data_ode.stan
}

transformed data {
#include stan_chunks/transformed_data.stan
}

parameters{
#include stan_chunks/parameters.stan
}

transformed parameters {
#include stan_chunks/transformed_parameters.stan
  y = integrate_ode_emp(init, t0, ts, theta, STEP_SIZE, INTERP_R, INTERP_A,
      x_r, x_i);
#include stan_chunks/likelihood_and_prior.stan
}

model {
#include stan_chunks/model.stan
}

generated quantities{
#include stan_chunks/generated_quantities.stan
}
