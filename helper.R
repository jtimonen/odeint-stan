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

