#######################################################
#####    
#####   Makefile for non-linear recovery project
#####    
#######################################################


###   Test the agent-based model against ODEs and Stochastic implementation

comp_ode: comp_ode_stoch_ABM.R SIR_nlr_ABM.R SIR_nlr_ODE.R SIR_nlr_stoch.R
	Rscript comp_ode_stoch_ABM.R


###   Generation interval analysis

analysis: SIR_nlr_ABM.R nlr_analysis.R prm.csv
	Rscript nlr_analysis.R
