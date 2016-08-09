
summary.gi <- function(gi.abm, CIwidth) {
	# Weekly Averages 
	gi.abm$day <- floor(gi.abm$t.round)
	gi.abm$week <- 1+ gi.abm$day %/% 7 
	gi.weekly.mean <- ddply(gi.abm,'week',summarize,
							gi = mean(gi_bck.mean),
							gi.lo = quantile(x=gi_bck.mean, probs = 0.05),
							gi.hi = quantile(x=gi_bck.mean, probs = 0.95))
	return(list(gi.weekly.mean = gi.weekly.mean,
				gi.abm = gi.abm))
}
