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
