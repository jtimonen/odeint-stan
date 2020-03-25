# Author: Juho Timonen
require(ggplot2)
require(bayesplot)

# Compute number of needed steps
compute_n_steps <- function(t_data, step_size, t0){
  dT <- max(t_data)-t0
  n_steps <- ceiling(dT/step_size) + 1
  return(n_steps)
} 

# For each measurement time point, finds the two points (or just the left one)
# between which it falls on an equispaced grid determined by step_size
compute_i_left <- function(t_data, step_size, t0){
  if(step_size<=0){stop('step_size must be positive')}
  N <- length(t_data)
  i_left <- rep(0, N)
  for(i in 1:N){
    idx <- 0
    while(t0 + idx*step_size < t_data[i]){ idx <- idx + 1}
    i_left[i] <- idx
    if(idx==0){stop('this should not occur')}
  }
  return(i_left)
}

# Format data for Stan
create_stan_data <- function(y0, t_data, Y_data, P, solver, step_size=NULL,
                             rel_tol=NULL, abs_tol=NULL, 
                             ABS_TOL=1e-10, REL_TOL=1e-10, max_steps=NULL,
                             n_steps_per_timepoint=NULL){
  N <- dim(Y_data)[1]
  D <- dim(Y_data)[2]
  methods <- c("rk45_boost", "euler", "euler_fixed",
               "rk4", "rk4_fixed")
  method <- which(methods==solver)-1
  t0 <- 0
  if(method %in% c(0,1,3)){
    if(!is.null(step_size)){
      stop(paste0("do not specify step size when solver=", solver))
    }
    step_size <- 1 # won't have any effect
  }
  n_steps <- compute_n_steps(t_data, step_size, t0)
  i_left <- compute_i_left(t_data, step_size, t0)
  if(is.null(rel_tol)){rel_tol <- 1e-6}
  if(is.null(abs_tol)){abs_tol <- 1e-6}
  if(is.null(max_steps)){max_steps <- 1e6}
  if(is.null(n_steps_per_timepoint)){n_steps_per_timepoint <- 1}
  out <- list(N=N, D=D, P=P, y0=y0, t0=t0, t=t_data, y=Y_data, method=method,
              step_size=step_size, n_steps=n_steps, i_left=i_left,
              n_steps_per_timepoint = n_steps_per_timepoint,
              rel_tol=rel_tol, abs_tol=abs_tol, max_steps=max_steps,
              ABS_TOL=ABS_TOL, REL_TOL=REL_TOL)
  return(out)
}

# Convenience wrapper for a function obtained by expose_stan_functions
odeint <- function(solver, y0, t0, t, theta, step_size=NULL,
                   n_steps_per_timepoint=NULL, rel_tol=NULL, 
                   abs_tol=NULL, max_steps=NULL){
  x_r <- x_i <- rep(0, 0)
  if(is.null(step_size)){
    if(is.null(n_steps_per_timepoint)){
      # rk45
      if(is.null(rel_tol)){rel_tol <- 1e-6}
      if(is.null(abs_tol)){abs_tol <- 1e-6}
      if(is.null(max_steps)){max_steps <- 1e6}
      Y <- solver(y0, t0, t, theta, x_r, x_i, rel_tol, abs_tol, max_steps) 
    }else{
      # euler or rk4
      Y <- solver(y0, t0, t, theta, x_r, x_i, n_steps_per_timepoint)
    }
  }else{
    # euler_fixed or rk4_fixed
    n_steps <- compute_n_steps(t, step_size, t0)
    i_left <- compute_i_left(t, step_size, t0)
    Y <- solver(y0, t0, t, theta, x_r, x_i, step_size, n_steps, i_left)
  }
  Y <- matrix(unlist(Y), length(t), length(y0), byrow = TRUE)
  return(Y)
}

# Takes subsamples of theta
subsample_theta <- function(fit, n_samples){
  ext <- rstan::extract(fit, pars=c('theta'), inc_warmup=TRUE, permuted=FALSE)
  samples <- rbind(ext[,1,], ext[,2,], ext[,3,])
  samples <- samples[sample(nrow(samples)),]
  samples <- samples[1:n_samples,]
  return(samples)
}

# Obtains solutions on a smoother time grid
solutions <- function(data, odeint_method, THETA, t, 
                      step_size=NULL, n_steps_per_timepoint=NULL,
                      rel_tol=NULL, abs_tol=NULL, max_steps=NULL){
  K  <- dim(THETA)[1]
  Y_hat <- array(0, c(K, length(t), length(data$y0)))
  for(k in 1:K){
    theta <- THETA[k, ]
    y_hat <- odeint(odeint_method, data$y0, 0, t, theta, 
                    step_size, n_steps_per_timepoint,
                    rel_tol, abs_tol, max_steps)
    Y_hat[k,,] <- y_hat
  }
  return(Y_hat)
}

# Plot data
plot_data <- function(t_data, Y_data){
  par(mfrow=c(1,2), oma = c(0, 0, 2, 0), mar=c(4,2.5,2.5,2))
  for(j in 1:2){
    plot(t_data, Y_data[,j], 'o', pch=16, xlab='t', ylab='', main=paste0('y',j))
  }
  mtext("Simulated data", outer = TRUE, cex = 1.5)
}

# Traceplot for chains
plot_chains <- function(fit, highlight=NULL, ncol=2, font_size=14){
  ext <- rstan::extract(fit, inc_warmup=TRUE, permuted=FALSE)
  n_warmup <- dim(ext)[1]/2
  if(is.null(highlight)){
    plt <- mcmc_trace(ext, n_warmup = n_warmup, facet_args = list(ncol=ncol))
  }else{
    plt <- mcmc_trace_highlight(ext, highlight = highlight, 
                                n_warmup = n_warmup, facet_args = list(ncol=ncol))
  }
  plt <- plt + theme_minimal() 
  plt <- plt + ggtitle('Trace plots of each param. (shaded part is warmup)')
  plt <- plt + theme(text = element_text(size=font_size))
  return(plt)
}

# Visualize solutions against data
plot_solutions <- function(t_data, Y_data, t_hat, Y_hat, alpha){
  par(mfrow=c(1,2), oma = c(0, 0, 2, 0), mar=c(4,2.5,2.5,2))
  K <- dim(Y_hat)[1]
  for(j in 1:2){
    plot(t_data, Y_data[,j], pch=16, xlab='t', ylab='', main=paste0('y',j))
    for(k in 1:K){
      lines(t_hat, Y_hat[k,,j], col = rgb(0.9,0.1,0.1,alpha=alpha))
    }
    points(t_data, Y_data[,j], pch=16, col='black')
  }
  mtext("Posterior solutions", outer = TRUE, cex = 1.5)
}

# Visualize solutions
plot_traj <- function(t_hat, Y_hat, names){
  y_out1 <- Y_hat[[1]]
  K <- length(Y_hat)
  N <- dim(y_out1)[1]
  D <- dim(y_out1)[2]
  x <- rep(0, K*N*D)
  for(k in 1:K){
    i1 <- 1 + (k-1)*N*D
    i2 <- k*N*D
    x[i1:i2] <- as.numeric(Y_hat[[k]])
  }
  t <- rep(t_hat, K*D)
  Method <- as.factor(rep(names, each=N*D))
  vars <- paste0('x', 1:D)
  Variable <- as.factor(rep(rep(vars, each=N),K))
  df <- data.frame(Method, Variable, t, x)
  plt <- ggplot(df, aes(x=t,y=x,group=Method, color=Method, lty=Method))
  plt <- plt + geom_line(lwd=1) + facet_wrap(.~ Variable)
  plt <- plt + theme(legend.position = 'bottom')
  plt <- plt + ggtitle('ODE solutions with theta = [1,1]')
  return(plt)
}

