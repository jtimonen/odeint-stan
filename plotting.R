# Plot data
plot_data <- function(t_data, Y_data){
  par(mfrow=c(1,2), oma = c(0, 0, 2, 0), mar=c(4,2.5,2.5,2))
  for(j in 1:2){
    plot(t_data, Y_data[,j], 'o', pch=16, xlab='t', ylab='', main=paste0('y',j))
  }
  mtext("Simulated data", outer = TRUE, cex = 1.5)
}

# Traceplot for chains
plot_chains <- function(fit, highlight=NULL, n_warmup=500, ncol=2){
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
  mtext("Posterior solutions using RK45", outer = TRUE, cex = 1.5)
}
