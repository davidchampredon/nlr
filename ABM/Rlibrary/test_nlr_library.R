##################################################################
######
######    MINIMAL TEST FOR 'nlr' LIBRARY
######
######
##################################################################

library(ggplot2)

library(nlr,lib.loc = "./lib")


params <- list(R0 = 2.123,
			  latent_mean = 5,
			  infectious_mean = 5,
			  nE = 0,
			  nI = 1,
			  Kfct = 'one', # one, sqrt, lin, exp, pow, inv
			  Kfct_prm = c(1,-5))

simulParams <- list(horizon = 299.9,
					popSize = 10000,
					seed = 1234,
					init_I1 = 5,
					calc_WIW_Re = FALSE,
					doExact = FALSE,
					timeStepTauLeap = 0.1)

res <- run_SEIR_nlc(params, simulParams)

t <- res[['ts_time']]
prev <- res[['ts_prevalence']]
S <- res[['ts_S']]
par(mfrow=c(2,1))
plot(t,prev, typ='s',log='y')
plot(t,S, typ='s')

# gi <- res[['gi_bck_val']]
# gi.t <- res[['gi_bck_times']]
# 
# gi.t <- gi.t[gi>0]
# gi <- gi[gi>0]
# plot(gi.t,gi)
# hist(gi)
