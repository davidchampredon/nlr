
source('SIR_nlr_ODE.R')
source('SIR_nlr_stoch.R')
source('SIR_nlr_ABM.R')
libpath <- "./ABM/Rlibrary/lib/"
library(nlr,lib.loc = libpath)

save.to.file <- F#TRUE

### ==== Shared parameters ====

horizon = 1
pop.size <- 2E4   # <-- recommended: 2E4

mu <- 1/1000
gamma <- 365/20
R0 <- 2.5
eps <- mu/(gamma+mu)
beta <- R0*(gamma+mu)
infect.init = 1E-3

K.fct <- 'one'  
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
n.MC <- 9        # Monte carlo iterations

# Initial values
I0 = pop.size * infect.init   # <--- this needs to be large for cvg stochastic -> deterministic
x0 <- c(S = pop.size-I0,
		I = I0,
		R = 0)


prm.sir <- list(mu=mu, beta=beta, gama=gamma, K.fct = 'one')
prm     <- list(mu=mu, beta=beta, gama=gamma, K.fct = K.fct, K.prm = K.prm)


### ==== ABM parameters ====

params <- list(R0 = R0,
			   latent_mean = 0.00001,
			   infectious_mean = 1/gamma*365,
			   nE = 0,
			   nI = 1,
			   Kfct = K.fct, # one, sqrt, lin, exp, pow, inv
			   Kfct_prm = K.prm)

simulParams <- list(horizon = horizon*365,
					popSize = pop.size,
					seed = 1234,
					init_I1 = I0,
					calc_WIW_Re = FALSE,
					doExact = FALSE,
					timeStepTauLeap = 0.1)


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


### ==== ABM simulations ====

sim.abm <- run_MC_ABM(params, simulParams, 
					  nMC = n.MC, 
					  nCPU = 3, 
					  libpath = libpath)
time.abm <- sim.abm[['ts']]$t / 365
prev.abm <- sim.abm[['ts']]$prev.mean


# ==== Plots ====

if(save.to.file) pdf('plot_comp_ode.pdf', width=6, height = 12)

par(mfrow=c(2,1))
lwd <- 6

plot(x= sim$time, 
	 y= sim$I * pop.size, 
	 ylim = range(sim$I * pop.size, 
	 			 sim.sir$I* pop.size,
	 			 prev.abm),
	 xlab = 'time', ylab = 'prevalence',
	 typ = 'l', las = 1,
	 lwd = lwd,
	 col = 'grey70',
	 main = 'SIR non-linear recovery\nComparison ODE vs Stochastic vs ABM')
lines(sim.sir$time, sim.sir$I* pop.size, 
	  col='grey', lty=2)

lines(sim.stoch$roundtime, sim.stoch$i.mean, col='red', lwd = lwd/3)
# lines(sim.stoch$roundtime, sim.stoch$i.q5, col='orange', lwd = 1)
# lines(sim.stoch$roundtime, sim.stoch$i.q95, col='orange', lwd = 1)

lines(x = time.abm, y = prev.abm, lwd = lwd/3, col = 'green3')


legend(x='topright', 
	   legend = c('ODE', 'Stoch', 'ABM' ,'ODE linear recov'),
	   col = c('grey70','red', 'green','grey'),
	   lwd = c(lwd, lwd/3, lwd/3, 1),
	   lty = c(1,1,1,2))

# log scale:
plot(x= sim$time, 
	 y= sim$I * pop.size, 
	 ylim = range(sim$I * pop.size, sim.sir$I* pop.size),
	 xlab = 'time', ylab = 'prevalence',
	 lwd = lwd,
	 typ = 'l', las = 1,
	 col = 'grey70',
	 log = 'y',
	 main = 'log scale')
lines(sim.sir$time, sim.sir$I* pop.size, 
	  col='grey', lty=2)
lines(sim.stoch$roundtime, sim.stoch$i.mean, col='red', lwd = lwd/3)
lines(x = time.abm, y = prev.abm, lwd = lwd/3, col = 'green3')


if(save.to.file) dev.off()
