//
//  individual.h
//  Gillespie_SEmInR
//
//  Created by David CHAMPREDON on 2015-02-25.
//

#ifndef __Gillespie_SEmInR__individual__
#define __Gillespie_SEmInR__individual__

#include <stdio.h>
#include "dcTools.h"

class individual
{
	unsigned long	_ID;
	unsigned long	_infectorID;
	
	unsigned int	_infectiousStatus;
	unsigned int	_maxInfectiousStatus;
	
	double			_timeInfectiousnessStart;
	double			_timeInfectiousnessEnd;
	double			_timeDiseaseAcquisition;
	vector<double>	_timeDiseaseTransmit;
	
	double			_GIbck;
	vector<double>	_GIfwd;
	
	double			_infectiousDuration;

	
public:

	individual(){};
	
	void create();
	void clean();
	
	
	// ====== SET FUNCTIONS =======
	
	void set_ID(unsigned long i) {_ID=i;}
	
	void set_infectiousStatus(unsigned int s)    {_infectiousStatus=s;}
	void set_maxInfectiousStatus(unsigned int s) {_maxInfectiousStatus=s;}
	
	void set_timeDiseaseAcquisition(double x)  {_timeDiseaseAcquisition = x;}
	void set_timeDiseaseTransmit(double x)     {_timeDiseaseTransmit.push_back(x);}
	void set_timeInfectiousnessStart(double x) {_timeInfectiousnessStart = x;}
	void set_timeInfectiousnessEnd(double x)   {_timeInfectiousnessEnd = x;}
	
	void set_infectorID(unsigned long i) {_infectorID=i;}
	
	void set_GIbck(double x)      {_GIbck = x;}
	void set_GIfwd_incr(double x) {_GIfwd.push_back(x);}
	
	void set_infectiousDuration(double x) {_infectiousDuration = x;}
	
	
	// ====== GET FUNCTIONS =======

	unsigned int		get_infectiousStatus()      {return _infectiousStatus;}
	unsigned long		get_ID()                    {return _ID;}
	unsigned long		get_infectorID()            {return _infectorID;}
	
	double				get_timeDiseaseAcquisition(){return _timeDiseaseAcquisition;}
	vector<double>		get_timeDiseaseTransmit()   {return _timeDiseaseTransmit;}
	
	double				get_timeInfectiousnessStart() {return _timeInfectiousnessStart;}
	double				get_timeInfectiousnessEnd()   {return _timeInfectiousnessEnd;}
	
	
	double				get_GIbck() {return _GIbck;}
	vector<double>		get_GIfwd() {return _GIfwd;}
	
	double				get_infectiousDuration()    {return _infectiousDuration;}
	
	
	// ====== EPIDEMIC =======
	
	void incrementStatus();
	
	
	// ====== HELPERS =======
	
	void displayInfo();
	
};




#endif /* defined(__Gillespie_SEmInR__individual__) */
