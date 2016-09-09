//
//  simulator.h
//  Gillespie_SEmInR
//
//  Created by David CHAMPREDON on 2015-02-25.
//

#ifndef __Gillespie_SEmInR__simulator__
#define __Gillespie_SEmInR__simulator__

#include <stdio.h>
#include "individual.h"

class simulator
{
	vector<individual>	_indiv;
	
	unsigned long		_popSize;
	
	// keeps track of individuals in S,
	// in any I[k] and in R
	unsigned long		_count_S;
	unsigned long		_count_I;
	unsigned long		_count_R;
    
    vector<unsigned long> _id_indiv_S;   // keep track of susceptible individuals
    vector<unsigned long> _id_indiv_I;   // keep track of infectious individuals
	
	// number of E and I compartments
	int					_nE;
	int					_nI;
	
	double				_beta;
	vector<double>		_sigma;
	vector<double>		_gamma;
	
	
	string				_Kfct; // Non linear recovery function
	vector<double>		_Kfct_prm; // parameters for _Kfct
	
	// events times
	vector<double>			_time;
	
	// time series
	vector<unsigned long>	_prevalence;
	vector<unsigned long>	_cumIncidence;
	vector<unsigned long>	_nS;	// keep track of the number of susceptible
	vector<unsigned long>	_nR;	// keep track of the number of recovered
	
	// Realized effective reproductive number
	// 1st column: Disease acquisition date for infector
	// 2nd column: number of secondary cases generated by that infector
	// Number of row = population size
	Matrix				_Reff;
	
	// Matrix of 'Who Infected Who' at different time points
	// element WIW(i,j)=0 if indiv#i has NOT infected indiv#j
	// element WIW(i,j)=g if indiv#i infected indiv#j with generation interval=g
	vector<Matrix>		_WIW;
	
	// times when _WIW was recorded
	// (to avoid calculating it at each step,
	//  which would be memory expensive)
	vector<double>		_WIW_times;
	
	
	
	// Simulation functions:
	void			initialize(unsigned long initInfectious);
	
	double			eventRate_infection();
	double			eventRate_latencyProgress(unsigned int i);
	double			eventRate_infectiousProgress(unsigned int i);
	double			eventRate_all_latencyProgress();
	double			eventRate_all_infectiousProgress();
	
	unsigned int			drawEventType();
	double					drawNextEventInterval();
	vector<unsigned long>	drawNumberEvents_tauLeap(double timestep);
	
	void			actionOnEvent(unsigned int eventType, double time_event);
	
	
	
public:
	
	simulator(){};
	simulator(double beta,vector<double>sigma,vector<double>gamma,
			  unsigned long popSize,
			  int nE, int nI);
	
	
	// ===== SET FUNCTIONS =====
	
	void			set_Kfct(string x)             {_Kfct = x;}
	void			set_Kfct_prm(vector<double> x) {_Kfct_prm = x;}
	
	void			set_GIbck(unsigned long IDindiv, double gi);
	void			set_GIfwd(unsigned long IDindiv, double gi);
	
	void			set_infectiousStatus(unsigned long IDindiv, unsigned int s);
	void			set_infectorID(unsigned long IDinfectee, unsigned long IDinfector);
	
	void			set_timeDiseaseAcquisition(unsigned long IDinfectee, double t);
	void			set_timeDiseaseTransmit(unsigned long IDinfector, double t);
	
	void			set_timeInfectiousnessStart(unsigned long IDinfector, double t);
	void			set_timeInfectiousnessEnd(unsigned long IDinfector, double t);
	void			set_infectiousDuration(unsigned long IDinfectious);
	
	void			calc_WIW(double t);
	void			calc_Reff();

	
	// ===== GET FUNCTIONS =====
	
	unsigned long			get_popSize(){return _popSize;}
	vector<unsigned long>	get_cumIncidence(){return _cumIncidence;}
	vector<unsigned long>	get_prevalence(){return _prevalence;}
	vector<unsigned long>	get_nS(){return _nS;}
	vector<unsigned long>	get_nR(){return _nR;}
	
	vector<double>			get_time(){return _time;}
	

	
	
	// ===== Simulation FUNCTIONS =====
	
	void run(double horizon,
			 unsigned long initInfectious,
			 bool calc_WIW_Re,
			 bool doExact,
			 double timeStep);

	void run__OLD(double horizon,
			 unsigned long initInfectious,
			 bool calc_WIW_Re);

	
	void run_exact(double horizon,
				   unsigned long initInfectious,
				   bool calc_WIW_Re);
	
	void run_tauLeap(double horizon,
					 double timestepSize,
					 unsigned long initInfectious,
					 bool calc_WIW_Re);
	
	void clean_start();
	
	
	
	// ===== OTHER FUNCTIONS =====
	
	vector< vector<double> > get_GIbck();
	vector< vector<double> > get_infectiousDuration();
	
	unsigned long	census_status(unsigned int a, unsigned int b); // counts individuals b/w infectious status a and b
	unsigned long	census_status(unsigned int a); // counts individuals of infectious status a
	unsigned long	census_S(); // counts susceptible individuals
	unsigned long	census_I(); // counts individuals in I[k] for all k
	
	bool			at_least_one_S_and_I();
	bool			at_least_one_S_and_E_or_I();
	
	vector<unsigned long> census(); // counts everyone in each status
	vector<unsigned long> census_not_R(); // counts everyone except 'R' (removed) individuals
	
	vector<unsigned long> census_ID(unsigned int a); // Retrieve IDs of all individuals of infectious status 'a'
	vector<unsigned long> census_ID_I(); // Retrieve IDs of all infectious individuals
	
	double	get_timeDiseaseAcquisition(unsigned long ID);
	
	
	// Save time series to file
	
	void save_prevalence(string filename);
	void save_cumIncidence(string filename);
	void save_nS(string filename);
	void save_nR(string filename);
	void save_Reff(string filename);
	void save_GIbck(string filename);
	void save_GIfwd(string filename);
	
	void displayInfo();
	
	
};


#endif /* defined(__Gillespie_SEmInR__simulator__) */
