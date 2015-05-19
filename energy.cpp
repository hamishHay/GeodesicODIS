#include "field.h"
#include "energy.h"
#include "mass.h"
#include "globals.h"
#include <math.h>

Energy::Energy(Mesh * mesh, int lat, int lon, Globals * Consts, Field * UVel, Field * VVel, Mass * MassField) : Field (mesh, lat, lon) 
{
	consts = Consts;
	u = UVel;
	v = VVel;
	mass = MassField;
};

void Energy::UpdateKinE(void) {
	//double cellMass = 0.;

	for (int i = 0; i < fieldLatLen; i++) {
		for (int j = 0; j < fieldLonLen; j++) {
			this->solution[i][j] = 0.5*mass->solution[i][j] * (pow(u->solution[i][j], 2) + pow(v->solution[i][j], 2));
		}
	}
};

void Energy::UpdateDtKinEAvg(void) {
	double kineticSum = 0; //Joules

	//Sum for poles
	//kineticSum += this->solution[0][0];
	//kineticSum += this->solution[this->fieldLatLen-1][0];
	
	for (int i = 0; i < fieldLatLen; i++) {
		for (int j = 0; j < fieldLonLen; j++) {
			kineticSum += this->solution[i][j];
		}
	}

	dtKinEAvg.push_back(kineticSum / (4 * pi*consts->radius.Value()*consts->radius.Value())); //Joules per meter^2


	//Automatically update dissipation - ensures vectors of samelength
	UpdateDtDissEAvg();
};

void Energy::UpdateDtDissEAvg(void) {
	dtDissEAvg.push_back(0);

	switch (dissMode) {
		//Linear dissipation
	case 0:
		dtDissEAvg[dtKinEAvg.size() - 1] = dtKinEAvg[dtKinEAvg.size() - 1] * consts->alpha.Value(); //Joules per meter
	}
};


void Energy::UpdateOrbitalKinEAvg(int inc) {	
	orbitKinEAvg.push_back(0);

	int pos = orbitKinEAvg.size()-1;
	
	for (int i = dtKinEAvg.size()-1; i > dtKinEAvg.size() - 1 - inc; i--) {
		orbitKinEAvg[pos] += dtKinEAvg[i];
	}

	orbitKinEAvg[pos] /= inc;//consts->period.Value();

	//Automatically update dissipation
	UpdateOrbitalDissEAvg();

};

void Energy::UpdateOrbitalDissEAvg(void) {
	orbitDissEAvg.push_back(0);

	switch (dissMode) {
		//Linear dissipation
	case 0:
		orbitDissEAvg[orbitKinEAvg.size() - 1] = orbitKinEAvg[orbitKinEAvg.size() - 1];// *consts->alpha.Value();
	}

};


