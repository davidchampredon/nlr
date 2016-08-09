###
### SIR WITH GILLESPIE ALGORITHM
###
### Author: David Champredon (2014-10-31)
###

library(ggplot2);theme_set(theme_bw())
library(plyr)
library(adaptivetau)  


##################################
#### STOCHASTIC GILLESPIE 
##################################

### STATE TRANSITIONS

transitions = list(c(S=1), # birth
				   c(S=-1, I=+1),  # infection
				   c(S=0, I=-1, R=+1),    # recovery
				   c(S=-1), # death S
				   c(I=-1), # death I
				   c(R=-1) # death R
)

### TRANSITIONS RATES

trans.rate <- function(x, params, t)
{
	birth <- params$mu
	infection <- params$beta*x["S"]*x["I"]*(x["S"]>=1)/(x["S"]+x["I"]+x["R"])
	
	if(params$K.fct == 'one')  k.fct <-  1
	if(params$K.fct == 'lin')  k.fct <-  x["I"] /(x["S"]+x["I"]+x["R"])
	if(params$K.fct == 'sqrt') k.fct <-  sqrt(x["I"]/(x["S"]+x["I"]+x["R"]))
	if(params$K.fct == 'exp')  k.fct <-  exp(-params$K.prm[1] * x["I"]/(x["S"]+x["I"]+x["R"]))
	if(params$K.fct == 'rev')  k.fct <-  params$K.prm[1] - x["I"]/(x["S"]+x["I"]+x["R"])
	if(params$K.fct == 'inv')  k.fct <-  1/(params$K.prm[1] + params$K.prm[2] * x["I"]/(x["S"]+x["I"]+x["R"]))
	if(params$K.fct == 'pow')  k.fct <-  ( x["I"] /(x["S"]+x["I"]+x["R"]) ) ^ (params$K.prm[1] -1)
	
	
	recovery <- params$gama * x["I"] * (x["I"]>=1) * k.fct 
	
	death.S <- params$mu * x["S"]*(x["S"]>=1)
	death.I <- params$mu * x["I"]*(x["I"]>=1)
	death.R <- params$mu * x["R"]*(x["R"]>=1)
	
	return( c(birth, infection, recovery,
			  death.S, death.I, death.R) )
}


run.stoch.mc <- function(do.adaptivetau, 
						 n.MC,
						 x0,
						 transitions,
						 trans.rate,
						 params.sir, 
						 params,
						 epsilon,
						 horizon) {
	
	# Data frame holding all MC iterations
	all.sim     <- data.frame(t=NULL,I=NULL,mc=NULL)
	all.sim.sir <- data.frame(t=NULL,I=NULL,mc=NULL)
	
	# Monte Carlo loop:
	for(mc in 1:n.MC)
	{
		if (!do.adaptivetau)
		{
			res.ATAU = ssa.exact(init.values = x0,
								 transitions, 
								 trans.rate, 
								 params = params.sir,
								 tf=horizon)
		}
		
		if(do.adaptivetau)
		{
			res.ATAU.sir = ssa.adaptivetau(init.values = x0,
										   transitions, 
										   trans.rate, 
										   params = params.sir,
										   tl.params = list(epsilon=epsilon),  
										   tf=horizon)
			
			res.ATAU = ssa.adaptivetau(init.values = x0,
									   transitions, 
									   trans.rate, 
									   params = params,
									   tl.params = list(epsilon=epsilon),  
									   tf=horizon)
		}
		
		tmp.sir = data.frame(t=res.ATAU.sir[,1], 
							 I=res.ATAU.sir[,"I"], 
							 mc=rep(mc,nrow(res.ATAU.sir)))
		
		tmp = data.frame(t=res.ATAU[,1], 
						 I=res.ATAU[,"I"], 
						 mc=rep(mc,nrow(res.ATAU)))
		
		all.sim.sir = rbind(all.sim.sir, tmp)
		all.sim = rbind(all.sim, tmp)
	}
	
	# Tidy simulations:
	
	grain = 3   # larger value = finer
	all.sim$roundtime = round(all.sim$t*grain,digits = 2)/grain
	
	sim.mean = ddply(all.sim,
					 .variables = c("roundtime"),
					 summarize,
					 i.mean = mean(I),
					 i.med = median(I),
					 i.q5 = quantile(I,0.05),
					 i.q95 = quantile(I,0.95),
					 i.q25 = quantile(I,0.25),
					 i.q75 = quantile(I,0.75))
	
	
	
	return(list(all.sim.sir = all.sim.sir, 
				all.sim     = all.sim,
				sim.mean    = sim.mean))
}


if(FALSE){
	do.adaptivetau = TRUE
	epsilon = 0.05  # default=0.05 ; larger=faster but less accurate
	
	set.seed(1234)
	n.MC = 10   # Monte carlo iterations
	horizon = 10
	
	# SIR parameters
	mu <- 1/80
	gama <- 365/20
	R0 <- 4
	beta <- R0*gama
	
	# Population size needs to be very large for 
	# convergence stochastic -> deterministic
	pop.size = 1E5
	
	# Initial values
	I0 = pop.size/1e3   # <--- this needs to be large for cvg stochastic -> deterministic
	x0 <- c(S=pop.size-I0,
			I=I0,
			R=0)
	
	
	
	# Data frame holding all MC iterations
	all.sim     <- data.frame(t=NULL,I=NULL,mc=NULL)
	all.sim.sir <- data.frame(t=NULL,I=NULL,mc=NULL)
	t0 <- proc.time()
	
	params.sir <- list(mu=mu, beta=beta, gama=gama, K.fct = 'one')
	params     <- list(mu=mu, beta=beta, gama=gama, K.fct = 'one', K.prm = c(0.0001 , 2))
	
	for(mc in 1:n.MC)
	{
		if (!do.adaptivetau)
		{
			res.ATAU = ssa.exact(init.values = x0,
								 transitions, 
								 trans.rate, 
								 params = list(beta=beta, gama=gama),
								 tf=horizon)
		}
		
		if(do.adaptivetau)
		{
			res.ATAU.sir = ssa.adaptivetau(init.values = x0,
										   transitions, 
										   trans.rate, 
										   params = params.sir,
										   tl.params = list(epsilon=epsilon),  
										   tf=horizon)
			
			res.ATAU = ssa.adaptivetau(init.values = x0,
									   transitions, 
									   trans.rate, 
									   params = params,
									   tl.params = list(epsilon=epsilon),  
									   tf=horizon)
		}
		
		tmp.sir = data.frame(t=res.ATAU.sir[,1], 
							 I=res.ATAU.sir[,"I"], 
							 mc=rep(mc,nrow(res.ATAU.sir)))
		
		tmp = data.frame(t=res.ATAU[,1], 
						 I=res.ATAU[,"I"], 
						 mc=rep(mc,nrow(res.ATAU)))
		
		all.sim.sir = rbind(all.sim.sir, tmp)
		all.sim = rbind(all.sim, tmp)
	}
	
	t.ATAU = proc.time()-t0
	
	### FIZZLES
	
	thres = 0.1
	sim.nofizz = ddply(all.sim,.variables = c("mc"),
					   summarize,
					   i.max = max(I))
	SIR.IMAX =   pop.size*0.01  #max(SIR.det$I)
	mc.nofizz = which(sim.nofizz$i.max>thres*SIR.IMAX)
	p.fizz = 1-length(mc.nofizz)/n.MC
	
	all.sim.nofizz = all.sim[all.sim$mc %in% mc.nofizz,]
	
	remove.fizzles <- FALSE
	if(remove.fizzles) all.sim = all.sim.nofizz
	
	
	##################################
	#### PLOT 
	##################################
	
	grain = 3   # larger value = finer
	all.sim$roundtime = round(all.sim$t*grain,digits = 2)/grain
	# all.sim.nofizz$roundtime = round(all.sim.nofizz$t*grain,digits = 2)/grain
	#unique(all.sim$roundtime)
	
	all.sim2 = all.sim[all.sim$t<horizon,]
	nrow(all.sim)-nrow(all.sim2)
	
	sim.mx = ddply(all.sim2,.variables = c("mc"),
				   summarize,
				   tmax = max(t))
	mc.noext = which(sim.mx$tmax>0.9*horizon)
	
	sim.mean = ddply(all.sim,.variables = c("roundtime"),
					 summarize,
					 i.mean = mean(I),
					 i.med = median(I),
					 i.q5 = quantile(I,0.05),
					 i.q95 = quantile(I,0.95),
					 i.q25 = quantile(I,0.25),
					 i.q75 = quantile(I,0.75))
	
	# sim.mean.nofizz = ddply(all.sim.nofizz,
	# 						.variables = c("roundtime"),
	# 						summarize,
	# 						i.mean = mean(I),
	# 						i.med = median(I),
	# 						i.q5 = quantile(I,0.05),
	# 						i.q95 = quantile(I,0.95),
	# 						i.q25 = quantile(I,0.25),
	# 						i.q75 = quantile(I,0.75))
	
	g = ggplot(all.sim)
	g=g+geom_point(aes(x=t, y=I),
				   alpha=0.10,
				   size=1.2) 
	g = g + ggtitle("SIR simulation")
	
	
	g = g + geom_ribbon(data=sim.mean, aes(x=roundtime,
										   y=i.mean,
										   ymin=i.q5,
										   ymax=i.q95),
						fill="orange",alpha=0.2)
	
	g = g + geom_ribbon(data=sim.mean, aes(x=roundtime,
										   y=i.mean,
										   ymin=i.q25,
										   ymax=i.q75),
						fill="orange",alpha=0.5)
	
	# g = g + geom_line(data=sim.mean.nofizz, aes(x=roundtime,y=i.mean),colour="green",size=2)
	g = g + geom_line(data=sim.mean, aes(x=roundtime,y=i.mean),colour="red",size=2)
	g = g + geom_line(data=sim.mean, aes(x=roundtime,y=i.med),colour="red")
	# g = g + geom_line(data=SIR.det, aes(x=time,y=I),colour="blue",size=2)
	plot(g)
	
	
	# Time performance for Gillespie adaptivetau
	print(t.ATAU[3])
	
}