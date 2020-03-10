# Creates a complete stan program code from three parts
# by adding a functions block to a base file
create_stan_code <- function(odeint_file, 
                             odefun_file, 
                             base_file = 'stan/base.txt',
                             n_rows_remove_base = 12){
  base <- readLines(base_file)
  base <- base[(n_rows_remove_base+1):length(base)] # remove first rows from base
  start <- 'functions{'
  fun1 <- readLines(odeint_file)
  fun2 <- readLines(odefun_file)
  end  <- '}'
  all  <- c(start, fun1, fun2, end, base)
  code <- paste(all, collapse="\n")
  return(code)
}

# Convenience wrapper for 'odeint' obtained by expose_stan_functions
odeint_wrap <- function(y0, t0, t, theta){
  x_r <- rep(0, 0)
  x_i <- rep(0, 0)
  Y   <- odeint(y0, t0, t, theta, x_r, x_i)
  N   <- length(t)
  Y   <- unlist(Y)
  Y   <- matrix(Y, N, length(Y)/N, byrow = TRUE)
  return(Y)
}

# Data simulation
simulate <- function(y0, t_data, theta, sigma){
  Y_data <- odeint_wrap(y0, 0, t_data, theta)
  N <- dim(Y_data)[1]
  D <- dim(Y_data)[2]
  Y_data <- Y_data + matrix(rnorm(n=N*D, sd=sigma), N, D)
  return(Y_data)
}

# Create data for Stan
create_stan_data <- function(y0, t_data, Y_data, P){
  N <- dim(Y_data)[1]
  D <- dim(Y_data)[2]
  t0 <- 0
  out <- list(N=N, D=D, P=P, y0=y0, t0=t0, t=t_data, y=Y_data)
  return(out)
}

