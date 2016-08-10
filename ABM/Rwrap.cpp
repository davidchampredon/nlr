//
//  Rwrap.cpp
//  NLR
//
//  Created by David CHAMPREDON on 2016-07-28.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//


#include <Rcpp.h>
#include "dcDataFrame.h"
#include "simulator.h"
#include "individual.h"
#include "mc.h"

using namespace Rcpp;



List dcDataFrameToRcppList(dcDataFrame df,
						   bool addRowName){
	/// Converts a dcDataFrame to a Rccp list
	/// (helper function)
	
	unsigned long ncol = df.get_colname().size();
	
	// Translate the 'dcDataFrame' into a R list
	// (convert to data frame in R, I don't know how to do it here in Rcpp):
	Rcpp::List rcpplist;
	// each column of the dcDataFrame is a list:
	for(int j=0; j<ncol; j++)
		rcpplist.push_back(df.get_value().extractColumn(j));
	// set the associated column names:
	if(!addRowName) rcpplist.attr("names") = df.get_colname();
	
	// insert row names (as an additional column,
	// don't know how to do differently)
	if(addRowName) {
		rcpplist.push_back(df.get_rowname());
		vector<string> cn = df.get_colname();
		cn.push_back("rowname");
		rcpplist.attr("names") = cn;
	}
	return rcpplist;
}


// [[Rcpp::export]]
List run_SEIR_nlc(List params, List simulParams){
	/// This runs a SEIR model implemented
	/// in this agent-based model.
	/// The goal is to compare its output to an ODE model.
	
	vector<unsigned int> res;
	
	
	// Model parameters:
	
	double R0				= params["R0"];
	double latent_mean		= params["latent_mean"];
	double infectious_mean	= params["infectious_mean"];
	int nE					= params["nE"];
	int nI					= params["nI"];
	string	Kfct			= params["Kfct"];
	vector<double> Kfct_prm	= params["Kfct_prm"];
	
	// Simulation parameters:
	
	double horizon			= simulParams["horizon"];
	unsigned long popSize	= simulParams["popSize"];
	unsigned long seed		= simulParams["seed"];
	
	
	unsigned long initInfectious	= simulParams["init_I1"];
	bool calc_WIW_Re				= simulParams["calc_WIW_Re"];
	bool doExact					= simulParams["doExact"];
	double timeStepTauLeap			= simulParams["timeStepTauLeap"];
	
	// parameter to get rid of (?)
	int jobnum = seed;
	string fileparam = "";
	
	
	// Derive other variables
	double sigma0	= 1.0 / latent_mean;
	double gamma0	= 1.0 / infectious_mean;
	double beta		= R0 * gamma0;
	
	
	vector<double> sigma(nE);
	vector<double> gamma(nI);
	for (int i=0; i<nE; i++) sigma[i]=sigma0 * nE;
	for (int i=0; i<nI; i++) gamma[i]=gamma0 * nI;
	
	// Simulation setup:
	
	simulator SIM(beta, sigma, gamma, popSize, nE, nI);
	
	SIM.set_Kfct(Kfct);
	SIM.set_Kfct_prm(Kfct_prm);
	

	// Call C++ function
	
	MC_run(SIM, 1, horizon,
		   initInfectious,
		   jobnum, fileparam, calc_WIW_Re,
		   doExact, timeStepTauLeap);

	SIM.displayInfo();
	
	// Retrieve all results from simulation:
	
	vector< vector<double> > gi = SIM.get_GIbck();
	vector<unsigned long> prev = SIM.get_prevalence();
	vector<unsigned long> nS = SIM.get_nS();
	vector<unsigned long> nR = SIM.get_nR();
	
	// populations:
//	dcDataFrame pop_final = sim.get_world()[0].export_dcDataFrame();
	// epidemic time series
//	dcDataFrame ts = sim.timeseries();

	
	// Return R-formatted result:
	return List::create(Named("gi_bck_times") = gi[0],
						Named("gi_bck_val") = gi[1],
						Named("ts_time") = SIM.get_time(),
						Named("ts_prevalence") = prev,
						Named("ts_S") = nS,
						Named("ts_R") = nR
						);
}

