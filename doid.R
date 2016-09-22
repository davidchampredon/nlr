source('SIR_nlr_ODE.R')
library("deSolve")


doid.ode <- function(t,x, parms){
	
	S <- x[1] 
	I <- x[2]
	U <- x[3]
	
	with(as.list(parms),{ 
		dS <- mu - beta*S*I - mu*S
		
		if(K.fct=='exp')    K_I <- exp(-K.prm[1]*I)
		if(K.fct=='affpow') K_I <- (K.prm[1] + K.prm[2] * I )^(K.prm[3])
		
		dI <- beta*S*I - mu*I - gamma * I * K_I
		
		if(t <= alpha) {
			dU = 0
		}
		else{
			dU <- - gamma * K_I * U
		}
		out <- c(dS, dI, dU) 
		list(out)
	})
}


# ==== Results =====

horizon  <- 300 * 1
timestep <- 0.09

infect.pop <- 0.01

mu    <- 1/(999*365)
gamma <- 1/4
R0    <- 2.0

eps  <- mu/(gamma+mu)
beta <- R0*(gamma+mu)

K.fct <- 'affpow'
K.prm <- c(0.1, 1, 1)


parms.sir <- list( mu=mu, gamma=gamma, beta=beta, K.fct = 'one' )
parms     <- list( mu=mu, gamma=gamma, beta=beta, K.fct = K.fct, K.prm = K.prm)

inits  <- c( S = (1-infect.pop), I =infect.pop  )
inits2 <- c( S = (1-infect.pop), I =infect.pop , U = 1.00 )

dt <- seq(timestep,horizon,timestep)

sim.sir <- as.data.frame(lsoda(inits,  dt, SIR_nlr, parms=parms.sir))



calc.doid <- function(alpha, inits, dt, parms, tt.max) {
	tmp <- c(alpha = alpha)
	sim <- as.data.frame(lsoda(inits2, dt, doid.ode,  
							   parms = append(parms,tmp)) )

	
	idx <- which(dt>= alpha)
	tt0 <- dt[idx] - dt[idx[1]]
	doid <- sim$U[idx]
	doid <- doid/sum(doid)
	
	zz <- exp(-tt0 * gamma)
	zz <- zz / sum(zz)

	doid.mean <- sum(tt0 * doid )
	zz.mean   <- sum(tt0 * zz)
	
	tt     <- tt0[tt0<tt.max]
	doid   <- doid[tt0<tt.max]
	zz     <- zz[tt0<tt.max]
	
	return(list(alpha = alpha, 
				sim   = sim, 
				doid  = doid, 
				doid.mean = doid.mean,
				doid.sir = zz,
				doid.sir.mean = zz.mean,
				tsa   = tt ))
}


# alpha.vec <- seq(1,250, by = 10)
alpha.vec <-2^seq(0,8,by=0.5)
res <- list()

for(i in seq_along(alpha.vec)){
	print(paste(round(i/length(alpha.vec)*100),'%'))
	res[[i]] <- calc.doid(alpha = alpha.vec[i], 
						  inits=inits2, dt=dt, parms=parms,
						  tt.max = 100)
}

sim <- res[[1]]$sim

### ===== PLOTS ====

par(mfrow=c(2,3))

### Phase plot

plot(x = sim$S,
	 y = sim$I,
	 main = paste0("K(I) = ",K.fct),
	 xlim = range(sim$S,sim.sir$S),
	 ylim = range(sim$I,sim.sir$I),
	 lwd = 2, col='red',
	 xlab="Susceptibles", ylab="Infected",
	 typ="l", log='xy')

lines(x = sim.sir$S,
	  y =  sim.sir$I,
	  lwd = 1, lty=2)
grid()


### Time series

plot(sim$time,sim$I, 
	 type="l", col="red", 
	 lwd=5,
	 main = paste0("K(I) = ",K.fct),
	 ylim = range(sim$I, sim.sir$I),
	 xlab="Time", ylab="Proportion of infectious",
	 cex=2)
lines(sim.sir$time,sim.sir$I, lty=2)
grid()

plot(sim$time,sim$I, 
	 type="l", col="red", 
	 lwd=5,
	 main = paste0("K(I) = ",K.fct),
	 ylim = range(sim$I, sim.sir$I),
	 xlab="Time", ylab="Proportion of infectious",
	 cex=2,
	 log = 'y')
lines(sim.sir$time,sim.sir$I, lty=2)
grid()


mx <- -99

doid.mean.plot <- numeric(length(res))
for(i in seq_along(res)){
	doid.mean.plot[i] <- res[[i]]$doid.mean
	mx <- max(mx,res[[i]]$doid.mean,res[[i]]$doid.sir.mean) 
}

for(i in seq_along(res)){
	
	n <- length(res)
	tt         <- res[[i]]$tsa
	doid       <- res[[i]]$doid
	doid.mean  <- res[[i]]$doid.mean
	alpha      <- res[[i]]$alpha
	zz         <- res[[i]]$doid.sir
	zz.mean    <- res[[i]]$doid.sir.mean
	
	col.sir <- 'gold'
	col.nlr <- rgb(0.2,0,1,i/n)
	if(i==1){
		plot(tt, doid, 
			 main = paste('DOI Distribution at',
			 			 res[[1]]$alpha,'-',res[[n]]$alpha,'days'),
			 typ='l', 
			 ylab = '', 
			 xlab = 'Time since disease acquisition',
			 xlim = range(tt,mx),
			 log = 'y',
			 las = 1,
			 yaxt = 'n',
			 col = col.nlr,
			 lwd=3)
		abline(v = doid.mean, col=col.nlr, lwd = 2, lty=2)
		lines(tt , zz, col=col.sir, lty=3)
		abline(v = zz.mean, col=col.sir, lty=3)
		grid()
	}
	if(i>1){
		lines(tt, doid, 
			  col = col.nlr,
			  lwd = 3)
		abline(v = doid.mean, col=col.nlr, lwd = 2, lty=2)
	}
	
}

plot(x = alpha.vec, 
	 y = doid.mean.plot,
	 main = 'Mean DOI',
	 pch = 16,
	 xlab = 'Time (alpha)',
	 ylab = 'mean DOI',
	 typ='o')
grid()


