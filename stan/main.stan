// Author: Juho Timonen
functions{
  
  // Edit this to be whatever system you want
  real[] odefun(real t, real[] y, real[] theta, data real[] x_r, data int[] x_i){
    real dydt[2];
    dydt[1] = 1*y[1] - theta[1]*y[1]*y[2];
    dydt[2] = theta[1]*y[1]*y[2] - theta[2] * y[2];
    return dydt;
  }
  
  // Vector version of odefun
  vector odefun_vec(real t, vector y, real[] theta, real[] x_r, int[] x_i){
    return to_vector(odefun(t, to_array_1d(y), theta, x_r, x_i));
  }
  
  // Linear interpolation from fixed time grid
  real[,] interpolate(real h, int[] i_left, real[] t_out, matrix y_grid){
    int N = num_elements(t_out);
    int D = rows(y_grid);
    matrix[D, N] y_out;
    for(i in 1:N){
      int idx = i_left[i];
      real d1 = t_out[i] - h*(idx-1);
      real d2 = h*(idx) - t_out[i];
      y_out[:,i] = y_grid[:,idx]*d2/(d1+d2) + y_grid[:,idx+1]*d1/(d1+d2);
    }
    return to_array_2d(transpose(y_out));
  }
  
  // Boost implementation of RK45 already in Stan
  real[,] rk45_boost(real[] y0, real t0, real[] t, real[] theta, real[] x_r, int[] x_i){
    int D = num_elements(y0);
    int N = num_elements(t);
    real y_hat[N, D];
    y_hat = integrate_ode_rk45(odefun, y0, t0, t, theta, x_r, x_i);
    return y_hat; 
  }

  // Forward Euler method
  real[,] euler(real[] y0, real t0, real[] t, real[] theta, real[] x_r, int[] x_i, real h, int n_steps, int[] i_left){
    int D = num_elements(y0);
    int N = num_elements(t);
    real y_hat[N, D];
    real t_curr = 0.0;
    vector[D] y_curr = to_vector(y0);
    matrix[D, n_steps] y_grid = rep_matrix(0.0, D, n_steps);
    y_grid[:,1] = y_curr;
    for(i in 2:n_steps){
      y_curr = y_curr + h*odefun_vec(t_curr, y_curr, theta, x_r, x_i);
      t_curr = t_curr + h;
      y_grid[:, i] = y_curr;
    }
    y_hat = interpolate(h, i_left, t, y_grid);
    return y_hat;
  }
  
    // A 4th order Runge-Kutta method
  real[,] rk4(real[] y0, real t0, real[] t, real[] theta, real[] x_r, int[] x_i, real h, int n_steps, int[] i_left){
    int D = num_elements(y0);
    int N = num_elements(t);
    real y_hat[N, D];
    real t_curr = 0.0;
    vector[D] y_curr = to_vector(y0);
    vector[D] k1;
    vector[D] k2;
    vector[D] k3;
    vector[D] k4;
    matrix[D, n_steps] y_grid = rep_matrix(0.0, D, n_steps);
    y_grid[:,1] = y_curr;
    for(i in 2:n_steps){
      k1 = h*odefun_vec(t_curr,       y_curr,        theta, x_r, x_i);
      k2 = h*odefun_vec(t_curr + h/2, y_curr + k1/2, theta, x_r, x_i);
      k3 = h*odefun_vec(t_curr + h/2, y_curr + k2/2, theta, x_r, x_i);
      k4 = h*odefun_vec(t_curr + h,   y_curr + k3,   theta, x_r, x_i);
      t_curr = t_curr + h;
      y_curr = y_curr + (k1 + 2*k2 + 2*k3 + k4)/6;
      y_grid[:, i] = y_curr;
    }
    y_hat = interpolate(h, i_left, t, y_grid);
    return y_hat;
  }
  
}

data{
  // Modeled data
  int<lower=1> N; // number of data points
  int<lower=1> D; // number of state variables
  int<lower=0> P; // number of ODE parameters (theta)
  real t[N];      // time measurements
  real y[N, D];   // state measurements
  real t0;        // initial time
  real y0[D];     // initial state
  
  // Solver properties
  int<lower=0,upper=2> method;
  real<lower=0> step_size;
  int<lower=1> n_steps;
  int<lower=0> i_left[N];
}

transformed data{
  real x_r[0];
  int x_i[0];
}

parameters{
  real<lower=0> theta[P];
  real<lower=0> sigma;
}

model{
  real y_hat[N, D];
  if(method==0){
    y_hat = rk45_boost(y0, t0, t, theta, x_r, x_i);
  }else if(method==1){
    y_hat = euler(y0, t0, t, theta, x_r, x_i, step_size, n_steps, i_left);
  }else if(method==2){
    y_hat = rk4(y0, t0, t, theta, x_r, x_i, step_size, n_steps, i_left);
  }else{
    reject("unknown method")
  }
  
  for(d in 1:D){
    target += normal_lpdf(y[:,d] | y_hat[:,d], sigma);
  }
  theta ~ normal(0, 5);
  sigma ~ normal(0, 5);
}
