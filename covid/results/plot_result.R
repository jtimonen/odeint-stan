
res <- readRDS('all_results.rds')
names <- c('euler', 'midpoint', 'adams', 'rk4')
colors <- c('firebrick', 'steelblue', 'gray20', 'green3')

# Get reference runtimes
t_rk45 <- mean(res$rk45$TIMES)/60
t_bdf <- mean(res$bdf$TIMES)/60


par(mfrow=c(3,1))
methods <- c("Euler", "Explicit midpoint", "2-step Adams-Bashforth", "4th order Runge-Kutta")

# Plot A
plot(c(-2,8), c(0.5,0.5), 'l', xlim=c(0, 6), ylim=c(-0.5, 4), col='orange',
     lty=2, xlab='Step size (days)', ylab='Pareto k')
grid()
for(j in 1:4){
  r <- res[[names[j]]]
  lines(r$h, r$K, col=colors[j])
  points(r$h, r$K, col=colors[j], pch=16)
}
legend(3.2, 4.1, lty=1, pch=16, legend=methods, col=colors, bty='n')

# Plot B
plot(c(-2,8), c(t_rk45, t_rk45), 'l', xlim=c(0, 6), ylim=c(-0.5, 105), col='red', lwd=1,
     xlab='Step size (days)', ylab='Average runtime (minutes)')
lines(c(-2,8), c(t_rk45, t_rk45), 'l', col = 'blue', lty=2, lwd=1)

grid()
for(j in 1:4){
  r <- res[[names[j]]]
  t <- rowMeans(r$TIMES)/60
  lines(r$h, t, col=colors[j])
  points(r$h, t, col=colors[j], pch=16)
}
legend(1.2, 70.8, lty=1, pch=16, legend=methods, col=colors, bty='n')
legend(3.0, 94.8, lty=c(1,2), legend=c("RK45 (Stan default tolerances)", "BDF (Stan default tolerances)"), 
       col=c("red", "blue"), lwd=1, bty='n')

# Plot C
plot(c(-2,120), c(0.5,0.5), 'l', xlim=c(0, 60), ylim=c(-0.5, 4), col='orange',
     lty=2, xlab='Average runtime (minutes)', ylab='Pareto k')
grid()
for(j in 1:4){
  r <- res[[names[j]]]
  t <- rowMeans(r$TIMES)/60
  lines(t, r$K, col=colors[j])
  points(t, r$K, col=colors[j], pch=16)
}
legend(33, 4.0, lty=1, pch=16, legend=methods, col=colors, bty='n')
