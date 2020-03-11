# Plot data
plot_data <- function(t_data, Y_data, connect=TRUE){
  N <- dim(Y_data)[1]
  D <- dim(Y_data)[2]
  Y_data <- as.numeric(Y_data)
  t_data <- rep(t_data, D)
  var <- as.factor(rep(paste("var", c(1:D)), each=N))
  DF <- data.frame(var, t_data, Y_data)
  colnames(DF) <- c("var", "t", "y")
  plt <- ggplot(DF, aes(x=t, y=y,group=var)) + facet_wrap(.~var)
  plt <- plt + geom_point() + theme_minimal()
  if(connect){
    plt <- plt + geom_line()
  }
  return(plt)
}

# A color selector
get_color <- function(idx){
  if(idx==1){
    color <- 'firebrick3'
  }else if(idx==2){
    color <- 'steelblue3'
  }else{
    stop("invalid idx")
  }
  return(color)
}

# Traceplot for chains
plot_chains <- function(fit, highlight=NULL, n_warmup=1000, ncol=2){
  ext <- rstan::extract(fit, inc_warmup=TRUE, permuted=FALSE)
  if(is.null(highlight)){
    plt <- mcmc_trace(ext, n_warmup = n_warmup, facet_args = list(ncol=ncol))
  }else{
    plt <- mcmc_trace_highlight(ext, highlight = highlight, 
                                n_warmup = n_warmup, facet_args = list(ncol=ncol))
  }
  plt <- plt + theme_minimal() 
  plt <- plt + ggtitle('Trace plots of each param. (shaded part is warmup)')
  return(plt)
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
