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
plot_chains <- function(fit, highlight=NULL){
  if(is.null(highlight)){
    plt <- mcmc_trace(fit)
  }else{
    plt <- mcmc_trace_highlight(fit, highlight = highlight)
  }
  return(plt)
}
