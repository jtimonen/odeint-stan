library(rstan)
library(tidyverse)
library(lubridate)
library(readxl)
source("model13/data_management_china.R")
source("model13/data_stan_model13_china.R") # creates data_list_model13
source("postproc.R")

model <- stan_model(file = 'model13/model13_bdf.stan')

fit <- sampling(object = model,
                data = data_list_model13,
                iter = 40,
                control = list(adapt_delta=0.8, max_treedepth=10),
                init = 0.5,
                chains = 1,
                cores = 1,
                refresh = 10,
                seed = 123)

# tmax = 42
S <- extract_compartment(fit, 'S') # 9383112
E <- extract_compartment(fit, 'E') # 115.3819
print(S[2,3,4])
print(E[2,3,4])
