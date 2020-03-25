# COMPARING LOG PROBS WITH DIFFERENT TOLERANCES

require(rstan)
source('functions.R')
rstan_options(auto_write = TRUE)
color_scheme_set("red")

model <- stan_model('stan/main.stan')
expose_stan_functions(model)


## 1. Simulating data
set.seed(27342)
t0     <- 0
y0     <- c(1,2)
t_data <- seq(0.5,12,by=0.5)
theta  <- c(1,1)     # true parameter values
sigma  <- 0.2        # noise std
N      <- length(t_data)
D      <- length(y0)
P      <- length(theta)

# Create and plot data
noise <- matrix(rnorm(n=N*D, sd=sigma), N, D)
Y_data <- odeint(rk45_boost, y0, t0, t_data, theta) + noise
plot_data(t_data, Y_data)

# Fit with different tolerances
stan_seed <- 2334
iter <- 2000
chains <- 4

TOL <- c(1e-6, 1e-4, 1e-2)
data1 <- create_stan_data(y0, t_data, Y_data, P, "rk45_boost",
                         abs_tol=TOL[1], rel_tol=TOL[1],
                         ABS_TOL=1e-10, REL_TOL=1e-10)

data2 <- create_stan_data(y0, t_data, Y_data, P, "rk45_boost",
                          abs_tol=TOL[2], rel_tol=TOL[2],
                          ABS_TOL=1e-10, REL_TOL=1e-10)

data3 <- create_stan_data(y0, t_data, Y_data, P, "rk45_boost",
                          abs_tol=TOL[3], rel_tol=TOL[3],
                          ABS_TOL=1e-10, REL_TOL=1e-10)

fit1  <- sampling(model, data1, seed = stan_seed, iter=iter, chains=chains)
fit2  <- sampling(model, data2, seed = stan_seed, iter=iter, chains=chains)
fit3  <- sampling(model, data3, seed = stan_seed, iter=iter, chains=chains)

#Extract
LP1 <- rstan::extract(fit1, pars=c("LP"))$LP
lp1 <- rstan::extract(fit1, pars=c("lp__"))$lp__
LP2 <- rstan::extract(fit2, pars=c("LP"))$LP
lp2 <- rstan::extract(fit2, pars=c("lp__"))$lp__
LP3 <- rstan::extract(fit3, pars=c("LP"))$LP
lp3 <- rstan::extract(fit3, pars=c("lp__"))$lp__

par(mfrow=c(2,2))
plot(exp(LP1-lp1))
plot(exp(LP2-lp2))
plot(exp(LP3-lp3))

# PSIS
log_ratios <- as.numeric(LP3-lp3)
IS <- psis(log_ratios)
