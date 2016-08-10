library(ggplot2)
library(nlr,lib.loc = "./ABM/Rlibrary/lib/")

source('SIR_nlr_ABM.R')
source('nlr_utils.R')


t0 <- as.numeric(Sys.time())

set.seed(1234)
save.to.file <- T

### ==== ABM parameters ====

horizon   <- 365
pop.size  <- 1E4
mu        <- 1/1000
gamma     <- 1/9

R0        <- 2.25

eps    <- mu/(gamma+mu)
beta   <- R0*(gamma+mu)

infect.init <- 1E-3

I0 <- pop.size * infect.init   

# K.fct <- 'one'   # one, sqrt, lin, exp, pow, inv

n.CPU <- 3
n.MC  <- 3 * n.CPU   # Monte carlo iterations


base.prm <- list(R0 = R0,
				 latent_mean = 3.123,
				 infectious_mean = 1/gamma,
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

wrap_abm <- function(X,
					 base.prm, 
					 simulParams, 
					 nMC = n.MC, 
					 nCPU = n.CPU) {
	
	params <- c(base.prm, X)
	sim.abm <- run_MC_ABM(params, 
						  simulParams, 
						  nMC = n.MC, 
						  nCPU = n.CPU)
	gi.abm <- sim.abm[['gi']]
	gi.abm <- gi.abm[gi.abm$gi_bck.mean>0,]
	z <- summary.gi(gi.abm, CIwidth =0.90)
	gi.weekly.mean <- z$gi.weekly.mean
	gi.abm <- z$gi.abm
	
	return(list(gi.abm = gi.abm,
				gi.weekly.mean = gi.weekly.mean))
}

res <- lapply(X           = K.list, 
			  FUN         = wrap_abm,
			  base.prm    = base.prm, 
			  simulParams = simulParams, 
			  nMC         = n.MC, 
			  nCPU        = n.CPU)

df.list <- list()
dfall.list <- list()

for (i in 1:length(res)) {
	tmp <- res[[i]]$gi.weekly.mean
	tmp$Kfct <- K.list[[i]]$Kfct
	df.list[[i]] <- tmp
	rm(tmp)
	tmp <- res[[i]]$gi.abm
	tmp$Kfct <- K.list[[i]]$Kfct
	dfall.list[[i]] <- tmp
}

df <- do.call('rbind',df.list)
dfall <- do.call('rbind',dfall.list)


### ==== PLOTS ====

if(save.to.file) pdf('plot_gi.pdf', width=12, height = 8)

g0 <- ggplot(df) 

g <- g0 + geom_line(aes(x=week, y=gi, colour = Kfct),
							size = 2)
g <- g + scale_y_log10()
g <- g + ggtitle('Backward GI Weekly mean')
plot(g)

g2 <- g + geom_ribbon(aes(x=week, ymin=gi.lo, ymax=gi.hi, fill = Kfct), 
					 alpha=0.2)
g2 <- g2 + facet_wrap(~Kfct)
plot(g2)


g.distrib <- ggplot(dfall) + geom_density(aes(x= gi_bck.mean, colour=factor(week)),
										  adjust = 2)
g.distrib <- g.distrib + facet_wrap(~Kfct)
g.distrib <- g.distrib + ggtitle("Temporal variation of backward GI distributions")
g.distrib <- g.distrib + scale_x_log10()
plot(g.distrib)


if(save.to.file) dev.off()

t1 <- as.numeric(Sys.time())

message(paste('Time elapsed:',round((t1-t0)/60, 2), 'min'))



