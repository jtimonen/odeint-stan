library(rstan)
library(tidyverse)
library(lubridate)
library(readxl)

# Get data
source("model13/data_management_china.R")
source("model13/data_stan_model13_china.R") # creates data_list_model13
source("postproc.R")

# Compile model
model <- stan_model(file = 'model13/model13_bdf.stan')

# Additional data for ode_integrate_bdf
data_list_model13$EPS      <- 1.0E-9
data_list_model13$abs_tol  <- 1.0E-10
data_list_model13$rel_tol  <- 1.0E-10
data_list_model13$max_iter <- 1.0E3

# Run sampling
fit <- sampling(object  = model,
                data    = data_list_model13,
                iter    = 40,
                control = list(adapt_delta=0.8, max_treedepth=10),
                init    = 0.5,
                chains  = 1,
                cores   = 1,
                refresh = 10)

