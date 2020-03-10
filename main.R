require(rstan)
source('helper.R')
c1 <- 'firebrick3'
c2 <- 'steelblue3'

# Specify system and integrator here
odeint_file <- 'stan/odefun_lotka-volterra.txt'
odefun_file <- 'stan/odeint_rk45-boost.txt'

# Create model code and compile model
code  <- create_stan_code(odeint_file, odefun_file)
model <- stan_model(model_code = code)

# Expose odefun() and odeint() to R
expose_stan_functions(model)

# Simulate
t_data <- seq(0.5,10,by=0.5)
y0     <- c(1,2)
theta  <- c(1,1,1)
Y_data <- odeint_wrap(y0, 0, t_data, theta)

# Plot
plot(t_data, Y_data[,1], 'o', lwd=2, pch=16, col = c1, xlab='t', ylab='y', main='Data')
lines(t_data, Y_data[,2], 'o', lwd=2, pch=16, col =  c2)
