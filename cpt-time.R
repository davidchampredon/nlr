x <- 100
pop <- c(10,  15, 20, 40, x)  # in thousands
dur <- c(1.1, 1.58,  2.22,  12.44, NA) # in minutes

dur.old <- c(1.8, 3.4,  6,  39, NA) # in minutes

estim.dur <- function(pop,dur,x) {
	L <- lm(log(dur)~pop)
	z <- L$coefficients
	durx <- exp(z[1]+z[2]*x)
	durx/60/24
}

durx <- estim.dur(pop,dur,x) 
durx.old <- estim.dur(pop,dur=dur.old,x) 
	


plot(pop,log(dur), typ='o',
	 ylim=c(z[1]/2, log(60*24*max(durx,durx.old))), 
	 lwd=2)
lines(pop,log(dur.old), col='lightgrey', lty=2, typ='o')

points(x,log(60*24*durx),col='red',pch=16,cex=2)
points(x,log(60*24*durx.old),col='orange',pch=16,cex=1)
grid()
# abline(a=z[1],b=z[2], col='red')


