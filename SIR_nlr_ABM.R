library(plyr)
library(snowfall)

run_MC_ABM <- function(params, simulParams, nMC, nCPU, libpath, timebucket=5) {
	
	sfInit(parallel = TRUE, cpu = nCPU)
	sfLibrary(nlr,lib.loc = libpath)
	
	snow.wrap <- function(i) {
		simulParams[['seed']] <- 123 + i
		sim.abm <- run_SEIR_nlr(params, simulParams)
		
		df.abm <- data.frame(t    = sim.abm[['ts_time']],
							 prev = sim.abm[['ts_prevalence']],
							 S    = sim.abm[['ts_S']],
							 mc   = i)
		
		df.gi <- data.frame(t      = sim.abm[['gi_bck_times']],
							gi_bck = sim.abm[['gi_bck_val']],
							mc     = i)
		
		df.infdur <- data.frame(t      = sim.abm[['infdur_times']],
								infdur = sim.abm[['infdur_val']],
								mc     = i)
		
		return(list(df.abm    = df.abm , 
					df.gi     = df.gi,
					df.infdur = df.infdur))
	}
	
	idx.apply <- 1:nMC
	
	sfExportAll()
	res <- sfSapply(idx.apply, snow.wrap, simplify = FALSE)
	sfStop()
	
	## Merge all MC iterations in one data frame:
	df.abm    <- list()
	df.gi     <- list()
	df.infdur <- list()
	
	for(i in 1:length(res)){
		df.abm[[i]]     <- res[[i]][[1]]
		df.gi[[i]]      <- res[[i]][[2]]
		df.infdur[[i]]  <- res[[i]][[3]]
	}
	
	all.ts     <- do.call('rbind',df.abm)
	all.gi     <- do.call('rbind',df.gi)
	all.infdur <- do.call('rbind',df.infdur)
	
	all.gi$t.round     <- round(all.gi$t     / timebucket, 2) * timebucket
	all.infdur$t.round <- round(all.infdur$t / timebucket, 2) * timebucket
	
	## Average across all monte carlo iterations:
	
	ts.s <- ddply(all.ts,'t',summarize,
				  prev.mean = mean(prev),
				  S.mean = mean(S))
	
	gi.s <- ddply(all.gi, 't.round', summarize,
				  gi_bck.mean = mean(gi_bck))
	
	infdur.s <- ddply(all.infdur, 't.round', summarize,
					  infdur.mean = mean(infdur))
	
	return( list(ts         = ts.s, 
				 gi         = gi.s,
				 infdur     = infdur.s,
				 ts.all     = all.ts,
				 gi.all     = all.gi,
				 infdur.raw = all.infdur))
}



