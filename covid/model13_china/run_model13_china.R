# needs data_list_model13
library(rstan)

#data_list_model13

model <- stan_model(file = 'models/model13.stan')

fit <- sampling(object = model,
                data = data_list_model13,
                warmup = 500,
                iter = 1000,
                control = list(adapt_delta=0.8, max_treedepth=10),
                init = 0.5,
                chains = 4,
                cores = 4,
                refresh = 20)

# tmax = 42