# Some currently unused functions

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


# Test timing
timing_test <- function(odeint_fun, t_data, y0, n_rep, h, fix_val=1){
  th <- seq(h,2,by=h)
  H  <- length(th)
  th1 <- rep(th, each=H)
  th2 <- rep(th, times=H)
  th3 <- rep(fix_val, H*H)
  THETA <- cbind(th1, th2, th3)
  runtimes <- time_solves(odeint, y0, t_data, THETA, n_rep=n_rep)
  df <- data.frame(th1, th2, runtimes)
  plt <- ggplot(df, aes(x=th1, y=th2, z=runtimes)) + geom_contour_filled()
  return(plt)
}
