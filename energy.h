#ifndef ENERGY_H
#define ENERGY_H

#include "vector"
#include "field.h"
#include "mesh.h"

class KineticEnergy : public Field {
private:
	//double orbitalGlobalAvg;
public:
	KineticEnergy(Mesh*,int,int);
	
	//vector of length n/(period/dt*1000), for average kinetic energy for n timesteps in one orbit.
	std::vector<double> timeStepKinAvg; 

	//vector of length m, for average kinetic energy at periapse for m orbits in simulation.
	std::vector<double> orbitKinAvg;

	void UpdateKinE(Field * u, Field * v);

	//Function finds globally averaged kinetic energy at the current timestep and appends it
	//to timeStepGlobalAvg
	void TimeStepKinEAvg(void); 

	//Function finds orbtially and globally averaged kinetic energy at end of the simulation for
	//the last whole orbit complete
	void OrbitalKinEAvg(int inc);

	//void FindGlobalDissAvg(void);

	//void FindOrbitalGlobalDissAvg(void);

};

#endif