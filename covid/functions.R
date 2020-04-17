# Create additional Stan data related to some ODE solver methods
add_interpolation_data <- function(data_list, h){
  t0 <- data_list$t0
  ts <- data_list$ts
  R  <- compute_R(t0, ts, h)
  A  <- compute_A(t0, ts, h, R)
  data_list$STEP_SIZE <- h
  data_list$INTERP_R  <- R
  data_list$INTERP_A  <- round(A, digits=14)
  return(data_list)
}

# Compute integers r_1, ..., r_n
compute_R <- function(t0, ts, h){
  n <- length(ts)
  R <- rep(0, n)
  for(i in 1:n){
    r <- 0
    while(t0 + r*h < ts[i]){
      r <- r + 1
    }
    R[i] <- r - 1
  }
  return(R)
}

# Compute multipliers a_1, ..., a_n
compute_A <- function(t0, ts, h, R){
  n <- length(ts)
  A <- rep(0, n)
  for(i in 1:n){
    D_i <- ts[i] - (t0 + R[i]*h)
    A[i] <- (h - D_i)/h
  }
  return(A)
}

# Get compartment
extract_compartment <- function(fit, name){
  comp <- paste0('comp_', name)
  X <- rstan::extract(fit, pars=comp)[[comp]]
  return(X)
}

# Plot compartment
plot_compartment <- function(fit, name='S', j=1, alpha = 0.2, color = 'steelblue4'){
  X  <- extract_compartment(fit, name)
  nS <- dim(X)[1]
  nT <- dim(X)[2]
  idx <- as.factor(rep(1:nS, each=nT))
  x <- rep(1:nT, nS)
  y <- as.numeric(t(X[,,j]))
  
  df <- data.frame(idx, x, y)
  p <- ggplot(df, aes(x=x, y=y, group=idx)) + geom_line(alpha = alpha, color = color)
  p <- p + ggtitle(paste0('Compartment ', name, ' (age group ', j, ')'))
  p <- p + xlab('Day') + ylab('Size')
  return(p)
}

# Helper function
get_samples <- function(param){
  samples <- as.vector(rstan::extract(fit, pars=param)[[param]])
  return(samples)
}
