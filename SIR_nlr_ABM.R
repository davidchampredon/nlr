library(plyr)
library(snowfall)

run_MC_ABM <- function(params, simulParams, nMC, nCPU, timebucket=5) {
	
	sfInit(parallel = TRUE, cpu = nCPU)
	sfLibrary(nlr,lib.loc = "./ABM/Rlibrary/lib/")
	
	snow.wrap <- function(i) {
		simulParams[['seed']] <- 123 + i
		sim.abm <- run_SEIR_nlc(params, simulParams)
		
		df.abm <- data.frame(t = sim.abm[['ts_time']],
								  prev = sim.abm[['ts_prevalence']],
								  S = sim.abm[['ts_S']],
								  mc = i)
		
		df.gi <- data.frame(t = sim.abm[['gi_bck_times']],
								 gi_bck = sim.abm[['gi_bck_val']],
								 mc = i)
		return(list(df.abm=df.abm , df.gi=df.gi))
	}
	
	idx.apply <- 1:nMC
	
	sfExportAll()
	res <- sfSapply(idx.apply, snow.wrap, simplify = FALSE)
	sfStop()
	
	## Merge all MC iterations in one data frame:
	df.abm <- list()
	df.gi <- list()
	for(i in 1:length(res)){
		df.abm[[i]] <- res[[i]][[1]]
		df.gi[[i]]  <- res[[i]][[2]]
	}
	
	all.ts <- do.call('rbind',df.abm)
	all.gi <- do.call('rbind',df.gi)
	all.gi$t.round <- round(all.gi$t / timebucket, 2) * timebucket
	
	
	## Average across all monte carlo iterations:
	
	ts.s <- ddply(all.ts,'t',summarize,
				  prev.mean = mean(prev),
				  S.mean = mean(S))
	gi.s <- ddply(all.gi,'t.round',summarize,
				  gi_bck.mean = mean(gi_bck))
	
	return( list(ts = ts.s, 
				 gi = gi.s,
				 ts.all = all.ts,
				 gi.all = all.gi))
}





run_MC_ABM_old <- function(params, simulParams, nMC) {
	
	#nMC <- 2 DEBUG
	
	df.abm <- list()
	df.gi <- list()
	
	for (i in 1:nMC) {
		simulParams[['seed']] <- 123 + i
		sim.abm <- run_SEIR_nlc(params, simulParams)
		
		df.abm[[i]] <- data.frame(t = sim.abm[['ts_time']],
								  prev = sim.abm[['ts_prevalence']],
								  S = sim.abm[['ts_S']],
								  mc = i)
		
		df.gi[[i]] <- data.frame(t = sim.abm[['gi_bck_times']],
								 gi_bck = sim.abm[['gi_bck_val']],
								 mc = i)
	}
	
	## Merge all MC iterations in one data frame:
	
	all.ts <- do.call('rbind',df.abm)
	all.gi <- do.call('rbind',df.gi)
	all.gi$t.round <- round(all.gi$t,1)
	
	
	## Average across all monte carlo iterations:
	
	ts.s <- ddply(all.ts,'t',summarize,
				  prev.mean = mean(prev),
				  S.mean = mean(S))
	gi.s <- ddply(all.gi,'t.round',summarize,
				  gi_bck.mean = mean(gi_bck))
	
	return( list(ts = ts.s, 
				 gi = gi.s,
				 ts.all = all.ts,
				 gi.all = all.gi))
}




