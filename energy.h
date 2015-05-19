#ifndef ENERGY_H
#define ENERGY_H

#include "vector"
#include "field.h"
#include "mesh.h"
#include "mass.h"
#include "globals.h"

/*
class IEnergy
{
public:
	virtual void Update() = 0;
};*/

class Energy : public Field {
private:
	Globals * consts;
	Field * u;
	Field * v;
	Mass * mass;

	
	int dissMode = 0; //dissipation type i.e., 0 for linear

	void UpdateDtDissEAvg(void);
	void UpdateOrbitalDissEAvg(void);

public:
	Energy(Mesh*, int, int, Globals *, Field *, Field *, Mass *);
	
	//vector of length n/(period/dt*1000), for average kinetic energy for n timesteps in one orbit.
	std::vector<double> dtKinEAvg;
	std::vector<double> dtDissEAvg;

	//vector of length m, for average kinetic energy at periapse for m orbits in simulation.
	std::vector<double> orbitKinEAvg;
	std::vector<double> orbitDissEAvg;

	void UpdateKinE(void);

	//Function finds globally averaged kinetic energy at the current timestep and appends it
	//to timeStepGlobalAvg
	void UpdateDtKinEAvg(void);
	

	//Function finds orbtially and globally averaged kinetic energy at end of the simulation for
	//the last whole orbit complete
	void UpdateOrbitalKinEAvg(int inc);
	
};
/*
class DissipatedEnergy : public Field {
public:
	DissipatedEnergy(Mesh *, int, int);


};*/

#endif