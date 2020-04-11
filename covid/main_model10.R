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





A <- extract_compartment(fit, 'R_GEN_')
B <- extract_compartment(fit, 'R')
diff <- as.vector(A - B)
print(max(diff^2))

var <- 'output_agedistr_deaths'
var_gen <- paste0(var, '_GEN_')
X <- rstan::extract(fit, pars = var)[[1]]
Y <- rstan::extract(fit, pars = var_gen)[[1]]
print(max(as.numeric(X - Y)^2))

