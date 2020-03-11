

# Convenience wrapper for a function obtained by expose_stan_functions
odeint_wrap <- function(odeint_method, y0, t0, t, theta){
  x_r <- x_i <- rep(0, 0)
  Y   <- odeint_method(y0, t0, t, theta, x_r, x_i)
  Y   <- matrix(unlist(Y), length(t), length(y0), byrow = TRUE)
  return(Y)
}

# Create data for Stan
create_stan_data <- function(y0, t_data, Y_data, P){
  N <- dim(Y_data)[1]
  D <- dim(Y_data)[2]
  out <- list(N=N, D=D, P=P, y0=y0, t0=0, t=t_data, y=Y_data)
  return(out)
}

# Time ode solving using different parameter values
time_solves <- function(odeint_fun, y0, t, THETA, n_rep=10){
  x_r <- rep(0, 0)
  x_i <- rep(0, 0)
  K <- dim(THETA)[1]
  t0 <- 0
  TIMES <- rep(0, K)
  for(k in 1:K){
    theta <- THETA[k,]
    start_time <- Sys.time()
    for(i in 1:n_rep){
      y_hat <- odeint_fun(y0, t0, t, theta, x_r, x_i)
    }
    TIMES[k] <- as.numeric(Sys.time() - start_time)/n_rep
  }
  return(TIMES)
}
