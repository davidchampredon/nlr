//
//  main_unit.cpp
//  Gillespie_SEmInR
//
//  Created by David CHAMPREDON on 2015-03-04.
//


#include <stdlib.h>
#include <iostream>
#include "simulator.h"
#include "individual.h"
#include "mc.h"

int main(int argc, const char * argv[]) {
	
	system("pwd");
	system("date");
	
	// For performance monitoring
	// - do not delete -
	timeval tim;
	gettimeofday(&tim, NULL);
	double t1=tim.tv_sec+(tim.tv_usec/1000000.0);
	// ------------------------------------------
	
	
	// Read the job number
    int jobnum = 0; //atoi(argv[1]);
	
	
	// main simulation parameters:
	
	double horizon			= 360;
    unsigned long popSize	= 15000;
	double R0				= 2.0;
	double latent_mean		= 2.1;
	double infectious_mean	= 3.1;
    int nE					= 1;
	int nI					= 1;
	unsigned long mc_iter	= 1;
	int njobs				= 1;
	
	unsigned long initInfectious	= 3;
	bool calc_WIW_Re				= false;
    bool doExact					= 0;
	double timeStepTauLeap			= 0.2;
	
	int mc_job = int(mc_iter/njobs);
	
	
	// Derive other variables
	double sigma0	= 1/latent_mean;
	double gamma0	= 1/infectious_mean;
	double beta		= R0*gamma0;
	
	
	vector<double> sigma(nE);
	vector<double> gamma(nI);
	for (int i=0; i<nE; i++) sigma[i] = sigma0*nE;
	for (int i=0; i<nI; i++) gamma[i] = gamma0*nI;
	
	// Simulation
	
	simulator SIM(beta, sigma, gamma, popSize, nE, nI);
	
    SIM.set_Kfct("sqrt");
    vector<double> kprm {1.0, 2.0};
    SIM.set_Kfct_prm(kprm);
    
	MC_run(SIM, mc_job, horizon,initInfectious,
		   jobnum,calc_WIW_Re,
		   doExact,timeStepTauLeap);
	
    SIM.displayInfo();
    
	cout<<endl<<"--- Job #"<<jobnum<<" finished!"<<endl;
	
	
//    vector< vector<double> > gi = SIM.get_GIbck();
//    displayVector(gi[1]);
//    
//    vector< vector<double> > infdur = SIM.get_infectiousDuration();
//    displayVector(infdur[0]);
//    displayVector(infdur[1]);
//    
    vector<unsigned long> prev = SIM.get_prevalence();
    displayVector(prev);
    
	// --------------------------------------------------------------
	// COMPUTER TIME MONITORING - do not delete!
	
	gettimeofday(&tim, NULL);
	double t2=tim.tv_sec+(tim.tv_usec/1000000.0);
	
	int minutes = (int)((t2-t1)/60.0);
	double sec = (t2-t1)-minutes*60.0;
	cout << endl << " - - - Computational time : ";
	cout << minutes<<" min "<<sec<<" sec" << endl;
	
	// --------------------------------------------------------------
	
	
	
	return 0;
}
