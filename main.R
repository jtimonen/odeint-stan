# Packages and options
require(rstan)
require(ggplot2)
require(bayesplot)
source('functions.R')
source('plotting.R')
rstan_options(auto_write = TRUE)
color_scheme_set("red")

# Specify system and integrator here
odeint_file <- 'stan/odefun_lotka-volterra.txt'
odefun_file <- 'stan/odeint_rk45-boost.txt'

# Create model code and compile model
code  <- create_stan_code(odeint_file, odefun_file)
model <- stan_model(model_code = code)

# Expose odefun() and odeint() to R
expose_stan_functions(model)

# Simulate
set.seed(123)
t_data <- seq(0.5,12,by=0.5)
y0     <- c(1,2)
theta  <- c(1,1,1)
sigma  <- 0.2
N      <- length(t_data)
D      <- length(y0)
P      <- length(theta)
Y_data <- odeint_wrap(y0, 0, t_data, theta)
Y_data <- Y_data + matrix(rnorm(n=N*D, sd=sigma), N, D)

# Plot data
plot_data(t_data, Y_data)

# Fit
stan_seed <- 123
data <- create_stan_data(y0, t_data, Y_data, P)
fit  <- sampling(model, data, seed = stan_seed)

# Diagnose
plot_chains(fit)
