
res <- readRDS('all_results.rds')
names <- c('euler', 'midpoint', 'adams', 'rk4')
colors <- c('firebrick', 'steelblue', 'gray20', 'green3')

par(mfrow=c(2,1))
plot(c(-2,8), c(0.5,0.5), 'l', xlim=c(0, 6), ylim=c(-0.5, 4), col='orange',
     lty=2, xlab='Step size (days)', ylab='Pareto k')
grid()
for(j in 1:4){
  r <- res[[names[j]]]
  lines(r$h, r$K, col=colors[j])
  points(r$h, r$K, col=colors[j], pch=16)
}
legend(3.2, 3.8, lty=1, pch=16, legend=names, col=colors)


plot(c(-2,120), c(0.5,0.5), 'l', xlim=c(0, 60), ylim=c(-0.5, 4), col='orange',
     lty=2, xlab='Average runtime (minutes)', ylab='Pareto k')
grid()
for(j in 1:4){
  r <- res[[names[j]]]
  t <- rowMeans(r$TIMES)/60
  lines(t, r$K, col=colors[j])
  points(t, r$K, col=colors[j], pch=16)
}
legend(40, 3.8, lty=1, pch=16, legend=names, col=colors)