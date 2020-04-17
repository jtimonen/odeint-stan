functions {
#include stan_functions/switch.stan
#include stan_functions/SEIR.stan
#include stan_functions/interpolate.stan
#include stan_functions/extract.stan
#include stan_functions/likelihood.stan
#include stan_functions/prior.stan

  // Adams-Bashforth 2-step method
  real[,] integrate_ode_ab2(real[] y0, real t0, real[] ts, real[] theta,
      data real STEP_SIZE, data int[] INTERP_R, data real[] INTERP_A, 
      data real[] x_r, data int[] x_i){
        
    int d = size(y0);
    int n = size(ts);
    real x[n, d];
    int R_n = INTERP_R[n];
    vector[d] y[R_n+2];
    real t = t0;
    
    vector[d] y_tmp;
    vector[d] f0;
    vector[d] f1;
    y[1] = to_vector(y0);
    
    // Compute y_1 using explicit midpoint method
    y_tmp = y[1] + 0.5*STEP_SIZE*odefun(t, y[1], theta, x_r, x_i);
    y[2] = y[1] + STEP_SIZE*odefun(t + 0.5*STEP_SIZE, y_tmp, theta, x_r, x_i);
    
    // Compute y_j for j >= 2 using the 2-step method
    for(i in 1:R_n){
      f0 = odefun(t, y[i], theta, x_r, x_i);
      f1 = odefun(t + STEP_SIZE, y[i+1], theta, x_r, x_i);
      y[i+2] = y[i+1] + STEP_SIZE*(1.5*f1 - 0.5*f0);
      t = t + STEP_SIZE;
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
  y = integrate_ode_ab2(init, t0, ts, theta, STEP_SIZE, INTERP_R, INTERP_A,
      x_r, x_i);
#include stan_chunks/likelihood_and_prior.stan
}

model {
#include stan_chunks/model.stan
}

generated quantities{
#include stan_chunks/generated_quantities.stan
}
