libpath <-  "./ABM/Rlibrary/lib/"

library(ggplot2)
library(nlr,lib.loc = libpath)

source('SIR_nlr_ABM.R')
source('nlr_utils.R')


t0 <- as.numeric(Sys.time())

set.seed(1234)
save.to.file <- T

### ==== ABM parameters ====

prm <- read.csv('prm.csv',header = FALSE)

getprm <- function(prm, name) {
	return(prm[prm[,1]==name,2])
}

horizon         <- getprm(prm,'horizon')  # in days
pop.size        <- getprm(prm,'pop.size')  #1E4
latent.mean     <- getprm(prm,'dol_mean')  # in days
infectious.mean <- getprm(prm,'doi_mean')  # in days
R0              <- getprm(prm,'R0')

mu        <- 1/1000
gamma     <- 1/infectious.mean
eps       <- mu/(gamma+mu)
beta      <- R0*(gamma+mu)

infect.init <- 1E-3
I0 <- pop.size * infect.init   

n.CPU <- getprm(prm,'ncpu')
n.MC  <- getprm(prm,'nmccpu') * n.CPU   # Monte carlo iterations


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

powvec <- c(0, 0.5, 1, 2, 3)
np <- length(powvec)

for(i in 1:np){
	K.list[[i]] <- list(Kfct = 'affpow',
						Kfct_prm = c(0,1,powvec[i]))
}

# K.list[[np+1]] <- list(Kfct = 'exp',
# 					Kfct_prm = c(0))



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
	
	# Time series:
	ts.abm <- sim.abm[['ts']]
	
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
	
	return(list(ts.abm = ts.abm,
				gi.abm = gi.abm,
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
df.ts.list <- list()
dfall.list <- list()
df.infdur.list <- list()
dfall.infdur.list <- list()
df.infdur.raw <- list()

KfctMakeID <- function(df, K.list) {
	df$Kfct <- K.list[[i]]$Kfct
	df$Kfct_prm <- paste(K.list[[i]]$Kfct_prm,collapse='_')
	df$KfctID <- paste(df$Kfct, df$Kfct_prm, sep='_')
	return(df)
}

for (i in 1:length(res)) {
	# Time series:
	tmp.ts <- res[[i]]$ts.abm
	# tmp.ts$Kfct <- K.list[[i]]$Kfct
	# tmp.ts$Kfct_prm <- paste(K.list[[i]]$Kfct_prm,collapse='_')
	df.ts.list[[i]] <- KfctMakeID(tmp.ts, K.list)
	rm(tmp.ts)
	
	# GI:
	tmp <- res[[i]]$gi.weekly.mean
	# tmp$Kfct <- K.list[[i]]$Kfct
	# tmp$Kfct_prm <- paste(K.list[[i]]$Kfct_prm,collapse='_')
	# df.list[[i]] <- tmp
	df.list[[i]] <- KfctMakeID(tmp, K.list)
	
	rm(tmp)
	tmp <- res[[i]]$gi.abm
	# tmp$Kfct <- K.list[[i]]$Kfct
	# dfall.list[[i]] <- tmp
	dfall.list[[i]] <- KfctMakeID(tmp, K.list)
	
	# Infectious duration:
	tmp2 <- res[[i]]$infdur.weekly.mean
	# tmp2$Kfct <- K.list[[i]]$Kfct
	# tmp2$Kfct_prm <- paste(K.list[[i]]$Kfct_prm,collapse='_')
	# df.infdur.list[[i]] <- tmp2
	df.infdur.list[[i]] <- KfctMakeID(tmp2, K.list)
	
	rm(tmp2)
	tmp2 <- res[[i]]$infdur.abm
	# tmp2$Kfct <- K.list[[i]]$Kfct
	# dfall.infdur.list[[i]] <- tmp2
	dfall.infdur.list[[i]] <- KfctMakeID(tmp2, K.list)

	zz <- KfctMakeID(tmp2, K.list)
	
	tmp3 <- res[[i]]$infdur.raw
	# tmp3$Kfct <- K.list[[i]]$Kfct
	# tmp3$Kfct_prm <- paste(K.list[[i]]$Kfct_prm,collapse='_')
	# df.infdur.raw[[i]] <- tmp3
	df.infdur.raw[[i]] <- KfctMakeID(tmp3, K.list)
	rm(tmp3)
}

df    <- do.call('rbind',df.list)
df.ts <- do.call('rbind',df.ts.list)
dfall <- do.call('rbind',dfall.list)

df.infdur    <- do.call('rbind',df.infdur.list)
dfall.infdur <- do.call('rbind',dfall.infdur.list)
df.infdur.raw <- do.call('rbind',df.infdur.raw)

df.infdur <- subset(df.infdur, week>0)
dfall.infdur <- subset(dfall.infdur, week>0)


### ==== PLOTS ====

theme_set(theme_bw())

# Time series:
if(save.to.file) pdf('plot_ts.pdf', width=12, height = 8)
plot_ts(df.ts)
if(save.to.file) dev.off()


# Generation intervals:

if(save.to.file) pdf('plot_gi.pdf', width=12, height = 8)
plot_mean(df = df, type = 'gi', popsize = pop.size,mc = n.MC)
plot_distribution_mean(dfall = dfall, type = 'gi', popsize = pop.size,mc = n.MC)
if(save.to.file) dev.off()

# Infectious duration:

if(save.to.file) pdf('plot_infdur.pdf', width=12, height = 8)
plot_mean(df = df.infdur, type = 'infdur', popsize = pop.size,mc = n.MC)
plot_distribution_mean(dfall = dfall.infdur, type = 'infdur', popsize = pop.size,mc = n.MC)
plot_distribution_raw(df.infdur.raw, xmax=60)
plot_distribution_raw_time(df=df.infdur.raw, xmax=60, time.bucket=10)
if(save.to.file) dev.off()


save.image(file = 'simul.RData')

t1 <- as.numeric(Sys.time())

message(paste('Time elapsed:',round((t1-t0)/60, 2), 'min'))




