libpath <-  "./ABM/Rlibrary/lib/"

library(ggplot2)
library(nlr,lib.loc = libpath)

source('SIR_nlr_ABM.R')
source('nlr_utils.R')


t0 <- as.numeric(Sys.time())

set.seed(1234)
save.to.file <- T

### ==== ABM parameters ====

horizon   <- 365  # in days
pop.size  <- 1E4

infectious.mean <- 4.5  # in days
latent.mean     <- 3.5  # in days

R0        <- 3.00

mu        <- 1/1000
gamma     <- 1/infectious.mean
eps    <- mu/(gamma+mu)
beta   <- R0*(gamma+mu)
infect.init <- 1E-3
I0 <- pop.size * infect.init   

n.CPU <- 3
n.MC  <- 5 * n.CPU   # Monte carlo iterations


base.prm <- list(R0 = R0,
				 latent_mean     = latent.mean,
				 infectious_mean = infectious.mean,
				 nE = 1,
				 nI = 1)

simulParams <- list(horizon = horizon,
					popSize = pop.size,
					seed = 1234,
					init_I1 = I0,
					calc_WIW_Re = FALSE,
					doExact = FALSE,
					timeStepTauLeap = 0.1)


### ==== ABM simulations ====

K.list <- list()

K.list[[1]] <- list(Kfct = 'one',
					Kfct_prm = c(0))

K.list[[2]] <- list(Kfct = 'sqrt',
					Kfct_prm = c(0))

K.list[[3]] <- list(Kfct = 'lin',
					Kfct_prm = c(0))

K.list[[4]] <- list(Kfct = 'exp',
					Kfct_prm = c(2))

K.list[[5]] <- list(Kfct = 'pow',
					Kfct_prm = c(3))

K.list[[6]] <- list(Kfct = 'inv',
					Kfct_prm = c(2, -5))


# Define wrap function for parallel execution:

wrap_abm <- function(X,
					 base.prm, 
					 simulParams, 
					 nMC = n.MC, 
					 nCPU = n.CPU) {
	# Run simulation:
	params <- c(base.prm, X)
	sim.abm <- run_MC_ABM(params, 
						  simulParams, 
						  nMC = n.MC, 
						  nCPU = n.CPU,
						  libpath = libpath)
	
	# Generation interval results:
	gi.abm <- sim.abm[['gi']]
	gi.abm <- gi.abm[gi.abm$gi_bck.mean>0,]
	z <- summary.gi(gi.abm, CIwidth =0.90)
	gi.weekly.mean <- z$gi.weekly.mean
	gi.abm <- z$gi.abm
	
	# Infectiousness duration results:
	infdur.abm <- sim.abm[['infdur']]
	infdur.abm <- infdur.abm[infdur.abm$infdur.mean>0,]
	z <- summary.infdur(infdur.abm, CIwidth =0.90)
	infdur.weekly.mean <- z$infdur.weekly.mean
	infdur.abm <- z$infdur.abm
	
	infdur.raw <- sim.abm[['infdur.raw']]
	
	return(list(gi.abm = gi.abm,
				gi.weekly.mean = gi.weekly.mean,
				infdur.abm = infdur.abm,
				infdur.weekly.mean = infdur.weekly.mean,
				infdur.raw = infdur.raw
				)
	)
}

res <- lapply(X           = K.list, 
			  FUN         = wrap_abm,
			  base.prm    = base.prm, 
			  simulParams = simulParams, 
			  nMC         = n.MC, 
			  nCPU        = n.CPU)

df.list <- list()
dfall.list <- list()
df.infdur.list <- list()
dfall.infdur.list <- list()
df.infdur.raw <- list()

for (i in 1:length(res)) {
	# GI:
	tmp <- res[[i]]$gi.weekly.mean
	tmp$Kfct <- K.list[[i]]$Kfct
	df.list[[i]] <- tmp
	rm(tmp)
	tmp <- res[[i]]$gi.abm
	tmp$Kfct <- K.list[[i]]$Kfct
	dfall.list[[i]] <- tmp
	
	# Infectious duration:
	tmp2 <- res[[i]]$infdur.weekly.mean
	tmp2$Kfct <- K.list[[i]]$Kfct
	df.infdur.list[[i]] <- tmp2
	rm(tmp2)
	tmp2 <- res[[i]]$infdur.abm
	tmp2$Kfct <- K.list[[i]]$Kfct
	dfall.infdur.list[[i]] <- tmp2

	tmp3 <- res[[i]]$infdur.raw
	tmp3$Kfct <- K.list[[i]]$Kfct
	df.infdur.raw[[i]] <- tmp3
	rm(tmp3)
}

df    <- do.call('rbind',df.list)
dfall <- do.call('rbind',dfall.list)

df.infdur    <- do.call('rbind',df.infdur.list)
dfall.infdur <- do.call('rbind',dfall.infdur.list)
df.infdur.raw <- do.call('rbind',df.infdur.raw)

df.infdur <- subset(df.infdur, week>0)
dfall.infdur <- subset(dfall.infdur, week>0)


### ==== PLOTS ====

if(save.to.file) pdf('plot_gi.pdf', width=12, height = 8)
plot_mean(df = df, type = 'gi')
plot_distribution_mean(dfall = dfall, type = 'gi')
if(save.to.file) dev.off()

if(save.to.file) pdf('plot_infdur.pdf', width=12, height = 8)
plot_mean(df = df.infdur, type = 'infdur')
plot_distribution_mean(dfall = dfall.infdur, type = 'infdur')
plot_distribution_raw(df.infdur.raw, xmax=60)
plot_distribution_raw_time(df.infdur.raw, xmax=60)
if(save.to.file) dev.off()


t1 <- as.numeric(Sys.time())

message(paste('Time elapsed:',round((t1-t0)/60, 2), 'min'))




