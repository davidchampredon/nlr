x <- 50
pop <- c(10,  15, 20, 40, x)  # in thousands
dur <- c(1.8, 3.4,  6,  39, NA) # in minutes

L <- lm(log(dur)~pop)
z <- (L$coefficients)

durx <- exp(z[1]+z[2]*x)
durx/60/24


plot(pop,log(dur), typ='o', lwd=2,
	 ylim=c(z[1]/2, log(durx)))
grid()
abline(a=z[1],b=z[2], col='red')


