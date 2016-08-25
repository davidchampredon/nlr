
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

summary.infdur <- function(infdur.abm, CIwidth) {
	# Weekly Averages 
	infdur.abm$day <- floor(infdur.abm$t.round)
	infdur.abm$week <- 1+ infdur.abm$day %/% 7 
	infdur.weekly.mean <- ddply(infdur.abm,'week',summarize,
							infdur = mean(infdur.mean),
							infdur.lo = quantile(x=infdur.mean, probs = 0.05),
							infdur.hi = quantile(x=infdur.mean, probs = 0.95))
	
	return(list(infdur.weekly.mean = infdur.weekly.mean,
				infdur.abm         = infdur.abm))
}


plot_mean <- function(df, type){
	g0 <- ggplot(df) 
	g <- g0 + geom_line(aes_string(x='week', 
								   y=type, 
								   colour='Kfct'),
						size = 2)
	g <- g + scale_y_log10()
	if(type=='gi') title <- 'Backward GI (weekly mean)'
	if(type=='infdur') title <- 'Infectious duration (weekly mean)'
	g <- g + ggtitle(title)
	plot(g)
}


plot_distribution_mean <- function(dfall, type) {
	
	if(type=='gi') {
		xlabel <- 'gi_bck.mean'
		title <- "Temporal variation of backward GI distributions"
	}
	if(type=='infdur') {
		xlabel <- 'infdur.mean'
		title <- "Temporal variation of infectious duration distributions"
	}
	
	dfall$week <- as.factor(dfall$week)
	
	g.distrib <- ggplot(dfall) + geom_density(aes_string(x= xlabel, 
														 colour='week'),
											  adjust = 2)
	g.distrib <- g.distrib + facet_wrap(~Kfct)
	g.distrib <- g.distrib + ggtitle(title)
	g.distrib <- g.distrib + scale_x_log10()
	plot(g.distrib)
}

plot_distribution_raw <- function(df,xmax=80) {
	
	x <- subset(df, infdur>0 & t>0)
	
	g <- ggplot(x)+geom_density(aes(x=infdur, colour=Kfct),adjust=1, size=1.5)
	g <- g + coord_cartesian(xlim=c(0,xmax))
	
	xx <- seq(0,xmax,by=0.1)
	yy <- dexp(xx,rate=1/infectious.mean)
	dfcheck <- data.frame(xx=xx,yy=yy)
	
	g <- g + geom_line(data = dfcheck, aes(x=xx,y=yy), colour='black', linetype = 2)
	g <- g + ggtitle('Infectious duration distributions \n (at all calendar times)')
	plot(g)
}

plot_distribution_raw_time <- function(df,xmax=80) {
	
	x <- subset(df, infdur>0 & t>0)
	
	x$day <- floor(x$t.round)
	x$week <- 1+ x$day %/% 7 
	
	g <- ggplot(x)+geom_density(aes(x=infdur, colour=factor(week)),adjust=1, size=0.5) 
	g <- g + facet_wrap(~Kfct, scales = 'free')
	g <- g + coord_cartesian(xlim=c(0,xmax),ylim=c(0,0.2))
	g <- g + ggtitle('Infectious duration distributions')
	plot(g)

}


