
col1 <- "firebrick"
col2 <- 'gray20'
par(col.lab = col2)
par(col.axis = col1)

h <- 0.5;
x <- seq(0, 4*h, by = h)
y <- c(1.6, 1.75,1.85,1.88,1.85)
plot(x,y,'o', ylim=c(1.5,2), pch=16, ylab='y', xlab='t', col=col1, lty=2, bty="n",
     xaxt = 'n', yaxt = 'n')
lines(x, rep(1.5, 5), col=col1, lty=2)
for(i in 1:5){
  lines(c(x[i], x[i]), c(1.5, y[i]), col=col1, lty=2)
}

x_dat <- c(0.2, 0.65, 0.85, 1.75)
y_dat <- approx(x, y, x_dat)$y
points(x_dat, y_dat, pch=15, col=col2)
for(i in 1:4){
  lines(c(x_dat[i], x_dat[i]), c(1.5, y_dat[i]), lty=1, col=col2)
}
#points(x_dat, rep(1.5, 4), pch=15, col=col2)

xlab <- c(expression(t[0]),
          expression(t[0] + h),
          expression(t[0] + 2*h),
          expression(t[0] + 3*h),
          expression(t[0] + 4*h))
axis(side=1,at=c(0,h,2*h,3*h, 4*h),labels=xlab, col=col1, tick=FALSE)

legend(0, 1.93, legend=c(expression(y[solution]), expression(y[output])),
       col=c(col1, col2), lty=c(2,1))
