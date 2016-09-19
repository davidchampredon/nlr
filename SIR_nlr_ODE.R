library("deSolve")


SIR_nlr <- function(t, x, parms)
{
	S <- x[1] 
	I <- x[2]
	
	with(as.list(parms),{ 
	dS <- mu - beta*S*I - mu*S
	
	if(K.fct=='one')    dI <- beta*S*I - mu*I - gamma * I
	if(K.fct=='sqrt')   dI <- beta*S*I - mu*I - gamma * I * sqrt(I)
	if(K.fct=='lin')    dI <- beta*S*I - mu*I - gamma * I * I
	if(K.fct=='rev')    dI <- beta*S*I - mu*I - gamma * I * (K.prm[1] - I)
	if(K.fct=='exp')    dI <- beta*S*I - mu*I - gamma * I * exp(-K.prm[1]*I)
	if(K.fct=='inv')    dI <- beta*S*I - mu*I - gamma * I / (K.prm[1] + K.prm[2]*I)
	if(K.fct=='pow')    dI <- beta*S*I - mu*I - gamma * I * I^(K.prm[1]-1)
	
	out <- c(dS, dI) 
	list(out)
})
} 




if(FALSE){
	
	horizon = 3
	timestep = 1/365
	
	infect.pop = 0.01
	
	mu <- 1/80
	gamma <- 365/20
	R0 <- 1.8
	eps <- mu/(gamma+mu)
	beta <- R0*(gamma+mu)
	
	K.fct <- 'sqrt'
	K.prm <- c(2,-5)
	
	parms.sir <- list( mu=mu, gamma=gamma, beta=beta, K.fct = 'one' )
	parms     <- list( mu=mu, gamma=gamma, beta=beta, K.fct = K.fct )
	
	inits <- c( S = (1-infect.pop), I =infect.pop  )
	dt <- seq(0,horizon,timestep)
	
	sim.sir <- as.data.frame(lsoda(inits, dt, SIR_nlr, parms=parms.sir))
	sim <- as.data.frame(lsoda(c(inits,1), dt, SIR_nlr, parms=parms))
	
	### CHECKS 
	
	# Simple SIR Endemic equilibrium
	S.ee <- 1/R0
	I.ee <- (1-1/R0)*eps
	
	
	### ===== PLOTS ====
	
	par(mfrow=c(1,3))
	
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
	points(x=S.ee,y=I.ee,col="green3",pch=16)
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
	
	plot(sim$time,sim$I, 
		 type="l", col="red", 
		 lwd=5,
		 main = paste0("K(I) = ",K.fct),
		 ylim = range(sim$I, sim.sir$I),
		 xlab="Time", ylab="Proportion of infectious",
		 cex=2,
		 log = 'y')
	lines(sim.sir$time,sim.sir$I, lty=2)
}