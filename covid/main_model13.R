library(rstan)
library(tidyverse)
library(lubridate)
library(readxl)
library(loo)

# Get data
source("model13/data_management_china.R")
source("model13/data_stan_model13_china.R") # creates data_list_model13
source("functions.R")

# Compile model
#model <- stan_model(file = 'model13/model13_bdf.stan')
#model <- stan_model(file = 'model13/model13_euler.stan')
#model <- stan_model(file = 'model13/model13_emp.stan')
#model <- stan_model(file = 'model13/model13_ab2.stan')
model <- stan_model(file = 'model13/model13_rk4.stan')

# Additional data for reference method ode_integrate_bdf
data_list_model13$EPS           <- 1.0E-9
data_list_model13$abs_tol_REF_  <- 1.0E-10
data_list_model13$rel_tol_REF_  <- 1.0E-10
data_list_model13$max_iter_REF_ <- 1.0E6

# Additional data for midpoint method
step_size <- 0.7
data_list_model13 <- add_interpolation_data(data_list_model13, step_size)

# Run sampling
fit <- sampling(object  = model,
                data    = data_list_model13,
                iter    = 1000,
                control = list(adapt_delta=0.8, max_treedepth=10),
                init    = 0.5,
                chains  = 1,
                cores   = 1,
                refresh = 10)

# PSIS-loo
lh1 <- get_samples('log_lik_na')
lh2 <- get_samples('log_lik_na_REF_')
pr1 <- get_samples('log_prior_na')
pr2 <- get_samples('log_prior_na_REF_')
post1 <- lh1 + pr1
post2 <- lh2 + pr2

out <- psis(post1 - post2)
print(out$diagnostics)
