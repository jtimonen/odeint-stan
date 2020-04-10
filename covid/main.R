library(rstan)
library(loo)
library(ggplot2)

dat   <- read_rdump('bb/covid.data.R')
model <- stan_model('full.stan')

fit <- sampling(object      = model,
                data        = dat, 
                iter        = 1000,
                control     = list(adapt_delta = 0.8, max_treedepth = 10),
                chains      = 3,
                save_warmup = TRUE)
#                init        = 0.5

PM <- FALSE
IW <- FALSE

prior <- as.vector(rstan::extract(fit, pars='prior', permuted = PM, inc_warmup = IW))
priorG <- as.vector(rstan::extract(fit, pars='prior_GEN_', permuted = PM, inc_warmup = IW))
lp__ <- as.vector(rstan::extract(fit, pars='lp__', permuted = PM, inc_warmup = IW))
lik <- as.vector(rstan::extract(fit, pars='lik', permuted = PM, inc_warmup = IW))
likG <- as.vector(rstan::extract(fit, pars='lik_GEN_', permuted = PM, inc_warmup = IW))

lp_GEN_ <- priorG + likG
#lik2 <- as.vector(rstan::extract(fit, pars='lik_2', permuted = PM, inc_warmup = IW))
#plot(lp__, prior + lik1)
#plot(lp__, prior + lik2)

plot(lp__, lp_GEN_)

out <- psis(as.vector(lp__ - lp_GEN_))

# Get compartments
S <- rstan::extract(fit, pars='comp_S')$comp_S
E <- rstan::extract(fit, pars='comp_E')$comp_E
I <- rstan::extract(fit, pars='comp_I')$comp_I
A <- rstan::extract(fit, pars='comp_A')$comp_A
R <- rstan::extract(fit, pars='comp_R')$comp_R
C <- rstan::extract(fit, pars='comp_C')$comp_C

extract_compartment <- function(fit, name){
  comp <- paste0('comp_', name)
  X <- rstan::extract(fit, pars=comp)[[comp]]
  return(X)
}

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

plot_compartment(fit, 'E', 3)



A <- extract_compartment(fit, 'R_GEN_')
B <- extract_compartment(fit, 'R')
diff <- as.vector(A - B)
print(max(diff^2))

var <- 'output_agedistr_deaths'
var_gen <- paste0(var, '_GEN_')
X <- rstan::extract(fit, pars = var)[[1]]
Y <- rstan::extract(fit, pars = var_gen)[[1]]
print(max(as.numeric(X - Y)^2))

