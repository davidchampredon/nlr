x <- 200
pop <- c(10,  15, 20, 40, 100, x)  # in thousands
dur <- c(1.1, 1.58,  2.22,  12.44, 120,NA) # in minutes for ONE MC per CPU (divide total time by # of MC per CPU)



estim.dur <- function(pop,dur,x) {
	L <- lm(log(dur)~pop)
	z <- L$coefficients
	durx <- exp(z[1]+z[2]*x)
	print(durx/60/24)
}

durx <- estim.dur(pop,dur,x) 

plot(pop,log(dur), typ='o',
	 ylim=c(0.1, log(60*24*max(durx,durx.old))), 
	 lwd=2)

points(x,log(60*24*durx),col='red',pch=16,cex=2)
grid()
# abline(a=z[1],b=z[2], col='red')


