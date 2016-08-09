
source('SIR_nlr_ODE.R')
source('SIR_nlr_stoch.R')

### ==== Common parameters ====

horizon = 1

pop.size <- 1E5

mu <- 1/80
gamma <- 365/20
R0 <- 1.8
eps <- mu/(gamma+mu)
beta <- R0*(gamma+mu)
infect.init = 1E-3

K.fct <- 'pow'
K.prm <- c(2,-5)


### ==== ODE parameters ====

timestep = 1/365
inits <- c( S = 1-infect.init, 
			I = infect.init  )
dt <- seq(0,horizon,timestep)


### ==== Stochastic parameters =====

set.seed(1234)
do.adaptivetau <- TRUE
epsilon <- 0.05  # default=0.05 ; larger=faster but less accurate
n.MC = 10   # Monte carlo iterations



# Initial values
I0 = pop.size * infect.init   # <--- this needs to be large for cvg stochastic -> deterministic
x0 <- c(S = pop.size-I0,
		I = I0,
		R = 0)




prm.sir <- list(mu=mu, beta=beta, gama=gamma, K.fct = 'one')
prm     <- list(mu=mu, beta=beta, gama=gamma, K.fct = K.fct, K.prm = K.prm)


### ==== ODE Simulations ====

sim.sir <- as.data.frame(lsoda(inits, dt, SIR_nlr, parms = prm.sir))
sim     <- as.data.frame(lsoda(inits, dt, SIR_nlr, parms = prm))
# Checks Simple SIR Endemic equilibrium
S.ee <- 1/R0
I.ee <- (1-1/R0)*eps


### ==== Stochastic Simulations ====

res.stoch <- run.stoch.mc(do.adaptivetau, 
						  n.MC,
						  x0,
						  transitions,
						  trans.rate,
						  params.sir = prm.sir, 
						  params = prm,
						  epsilon,
						  horizon)

sim.stoch <- res.stoch[['sim.mean']]


# ==== Plots ====

par(mfrow=c(2,1))
lwd <- 5

plot(x= sim$time, 
	 y= sim$I * pop.size, 
	 ylim = range(sim$I * pop.size, sim.sir$I* pop.size),
	 xlab = 'time', ylab = 'prevalence',
	 typ = 'l', las = 1,
	 lwd = lwd,
	 main = 'SIR non-linear recovery\nComparison ODE vs Stochastic')
lines(sim.sir$time, sim.sir$I* pop.size, 
	  col='grey', lty=2)

lines(sim.stoch$roundtime, sim.stoch$i.mean, col='red', lwd = lwd/2)
lines(sim.stoch$roundtime, sim.stoch$i.q5, col='orange', lwd = 1)
lines(sim.stoch$roundtime, sim.stoch$i.q95, col='orange', lwd = 1)

legend(x='topright', 
	   legend = c('ODE', 'Stoch', 'ODE linear recov'),
	   col = c('black','red','grey'),
	   lwd = c(lwd,lwd/2,1),
	   lty = c(1,1,2))

# log scale:
plot(x= sim$time, 
	 y= sim$I * pop.size, 
	 ylim = range(sim$I * pop.size, sim.sir$I* pop.size),
	 xlab = 'time', ylab = 'prevalence',
	 typ = 'l', las = 1,
	 lwd = lwd,
	 log = 'y',
	 main = 'log scale')
lines(sim.sir$time, sim.sir$I* pop.size, 
	  col='grey', lty=2)
lines(sim.stoch$roundtime, sim.stoch$i.mean, col='red', lwd = lwd/2)