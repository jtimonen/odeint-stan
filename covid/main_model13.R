library(rstan)
library(tidyverse)
library(lubridate)
library(readxl)
library(loo)

# Get data
source("model13/data_management_china.R")
source("model13/data_stan_model13_china.R") # creates data_list_model13
source("postproc.R")

# Compile model
#model <- stan_model(file = 'model13/model13_bdf.stan')
model <- stan_model(file = 'model13/model13_midpoint.stan')

# Additional data for ode_integrate_bdf
data_list_model13$EPS      <- 1.0E-9
data_list_model13$abs_tol  <- 1.0E-10
data_list_model13$rel_tol  <- 1.0E-10
data_list_model13$max_iter <- 1.0E6

# Additional data for midpoint method
data_list_model13$step_size <- 1.0

# Run sampling
fit <- sampling(object  = model,
                data    = data_list_model13,
                iter    = 1000,
                control = list(adapt_delta=0.8, max_treedepth=10),
                init    = 0.5,
                chains  = 1,
                cores   = 1,
                refresh = 10)

# Helper function
get_samples <- function(param){
  samples <- as.vector(rstan::extract(fit, pars=param)[[param]])
  return(samples)
}

# PSIS-loo
lh1 <- get_samples('log_lik_na')
lh2 <- get_samples('log_lik_na_GEN_')
pr1 <- get_samples('log_prior_na')
pr2 <- get_samples('log_prior_na_GEN_')
post1 <- lh1 + pr1
post2 <- lh2 + pr2

out <- psis(post1 - post2)
print(out$diagnostics)
