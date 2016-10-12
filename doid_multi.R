source('SIR_nlr_ODE.R')
library("deSolve")


doid.ode <- function(t,x, parms){
	# Implement ODEs of the cohort model.
	
	S <- x[1] 
	I <- x[2]
	U <- x[3]
	
	with(as.list(parms),{ 
		dS <- mu - beta*S*I - mu*S
		
		if(K.fct=='exp')    K_I <- exp(-K.prm[1]*I)
		if(K.fct=='affpow') K_I <- (K.prm[1] + K.prm[2] * I )^(K.prm[3])
		
		dI <- beta*S*I - mu*I - gamma * I * K_I
		
		# '1-U' is the _cumulative_ distribution
		# of the infectious period:
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



calc.doid <- function(alpha, inits, dt, parms, tt.max) {
	# Solve for the DOI distribution
	# using ODEs from the cohort model.
	
	tmp   <- c(alpha = alpha)
	gamma <- parms[['gamma']]
	sim   <- as.data.frame(lsoda(inits, dt, doid.ode,  
								 parms = append(parms,tmp)) )
	
	idx <- which(dt>= alpha)
	tt0 <- dt[idx] - dt[idx[1]]
	# (1-U) is the cumulative distribution of 
	# the infectious period, so the _density_
	# is its derivative (negative):
	U <- sim$U[idx]
	# derivative:
	deltat <- dt[2]-dt[1]
	doid   <- - diff(U)/deltat
	# normalize:
	doid <- doid/sum(doid)
	
	# check with simple SIR:	
	doid.sir <- exp(-tt0 * gamma)
	doid.sir <- doid.sir / sum(doid.sir)
	
	tt0.m <- tt0[1:length(doid)]
	doid.mean     <- sum(tt0.m * doid )
	doid.sir.mean <- sum(tt0 * doid.sir)
	
	tt        <- tt0[tt0<tt.max]
	doid      <- doid[tt0<tt.max]
	doid.sir  <- doid.sir[tt0<tt.max]
	
	return(list(alpha = alpha, 
				sim   = sim, 
				doid  = doid, 
				doid.mean = doid.mean,
				doid.sir = doid.sir,
				doid.sir.mean = doid.sir.mean,
				tsa   = tt ))
}





simul.one <- function(classic.prm, K.fct, K.prm) {
	
	# Parameters set-up:
	
	horizon <- classic.prm$horizon   
	timestep <- classic.prm$timestep 
	infect.pop <- classic.prm$infect.pop
	mu <- classic.prm$mu          
	gamma <- classic.prm$gamma    
	R0 <- classic.prm$R0          
	
	eps  <- mu/(gamma+mu)
	beta <- R0*(gamma+mu)
	dt <- seq(timestep,horizon,timestep)
	
	parms.sir <- list( mu=mu, gamma=gamma, beta=beta, K.fct = 'one' )
	parms     <- list( mu=mu, gamma=gamma, beta=beta, K.fct = K.fct, K.prm = K.prm)
	inits  <- c( S = (1-infect.pop), I =infect.pop  )
	inits2 <- c( S = (1-infect.pop), I =infect.pop , U = 1.00 )
	
	# solve simple SIR ODEs:
	sim.sir <- as.data.frame(lsoda(inits,  dt, SIR_nlr, parms=parms.sir))
	
	# solve cohort model for DOI distribution
	# at several calendar date ('alpha')
	
	alpha.vec <-2^seq(0,8,by=0.5)
	
	res <- list()
	for(i in seq_along(alpha.vec)){
		print(paste(round(i/length(alpha.vec)*100),'%'))
		res[[i]] <- calc.doid(alpha = alpha.vec[i], 
							  inits=inits2, dt=dt, parms=parms,
							  tt.max = 100)
	}
	
	return(list(res=res, 
				sim.sir=sim.sir, 
				alpha.vec=alpha.vec))
}


### Time series

plot.ts <- function(sim, sim.sir, K.fct, K.prm, log.plot='y'){
	u <- paste(c('a','b','c'),K.prm,sep='=',collapse = ' ; ')
	
	plot(sim$time,sim$I, 
		 type="l", col="red", 
		 lwd=5,
		 main = paste0("K = ",K.fct,'\n',u),
		 ylim = pmax(1e-9,range(sim$I, sim.sir$I)),
		 xlab="Time", ylab="Proportion of infectious",
		 cex=2,
		 log = log.plot)
	lines(sim.sir$time, sim.sir$I, lty=2, col='grey')
	grid()
}



pot.doid.alpha <- function(res,log.plot='y'){
	
	mx <- -99
	for(i in seq_along(res)){
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
		col.nlr <- rgb(0, 0, 1-i/n, i/n)
		if(i==1){
			plot(tt, doid, 
				 main = paste('DOI Distribution (log-scale) at \n day',
				 			 res[[1]]$alpha,'(lightest) to day',res[[n]]$alpha,'(darkest)'),
				 typ='l', 
				 ylab = '', 
				 xlab = 'Time since disease acquisition',
				 xlim = range(tt,mx),
				 log = log.plot,
				 las = 1,
				 #yaxt = 'n',
				 col = col.nlr,
				 lwd =2)
			# abline(v = doid.mean, col=col.nlr, lwd = 2, lty=2)
			lines(tt , zz, col=col.sir, lty=3)
			# abline(v = zz.mean, col=col.sir, lty=3)
			grid()
		}
		if(i>1){
			lines(tt, doid, 
				  col = col.nlr,
				  lwd = 3)
			# abline(v = doid.mean, col=col.nlr, lwd = 2, lty=2)
		}
	}
}


plot.doid.mean <-function(res,alpha.vec){
	doid.mean.plot <- numeric(length(res))
	for(i in seq_along(res)){
		doid.mean.plot[i] <- res[[i]]$doid.mean
	}
	
	plot(x = alpha.vec, 
		 y = doid.mean.plot,
		 main = 'Mean DOI',
		 pch = 16,
		 xlab = 'Time (alpha)',
		 ylab = 'mean DOI',
		 typ='o')
	grid()
}


# ==== Parameters =====

classic.prm <- list()
classic.prm[['horizon']]     <- 300 * 1
classic.prm[['timestep']]    <- 0.09
classic.prm[['infect.pop']]  <- 1e-5
classic.prm[['mu']]          <- 1/(999*365)
classic.prm[['gamma']]       <- 1/4
classic.prm[['R0']]          <- 2.0


K.fct  <- 'affpow'
powvec <- c(1, 0.5, 1, 2, 1, 1,-1)
a      <- c(1, 0,   0, 0, 1, 1, 1)
b      <- c(0, 1,   1, 1, 1,-1, 2)

K <-list()
for(q in seq_along(a)){
	K[[q]] <- list(K.fct, c(a[q], b[q], powvec[q]))
}


### ===== Results & Plots ====

t0<- as.numeric(Sys.time())

save.plot.file <- TRUE
if(save.plot.file) pdf('plot_doid.pdf', height = 20, width = 15)

par(mfrow=c(length(K),3))
Z <- list()
for(i in seq_along(K)){
	print(paste(i,'/',length(K)))
	K.fct <- K[[i]][[1]]
	K.prm <- K[[i]][[2]]
	Z[[i]] <- simul.one(classic.prm, K.fct, K.prm)	
	res       <- Z[[i]]$res
	sim       <- res[[1]]$sim
	sim.sir   <- Z[[i]]$sim.sir
	alpha.vec <- Z[[i]]$alpha.vec
	
	plot.ts(sim, sim.sir, K.fct, K.prm, log.plot='')
	pot.doid.alpha(res, log.plot='')
	plot.doid.mean(res,alpha.vec)
}
if(save.plot.file) dev.off()

t1<- as.numeric(Sys.time())

print(paste('Time elapsed:',round((t1-t0)/60,1),'min'))




