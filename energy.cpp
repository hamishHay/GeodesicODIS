#include "field.h"
#include "energy.h"
#include "mass.h"
#include "globals.h"
#include "outFiles.h"
#include <math.h>

#include <iostream>
#include <fstream>
#include <iomanip>

Energy::Energy(Mesh * mesh, int lat, int lon, Globals * Consts, Field * UVel, Field * VVel, Mass * MassField) : Field (mesh, lat, lon) 
{
	consts = Consts;
	u = UVel;
	v = VVel;
	mass = MassField;
};

void Energy::UpdateKinE(void) {
	//double cellMass = 0.;

	switch (consts->fric_type) {
		//Linear dissipation
	case LINEAR:
		for (int i = 0; i < fieldLatLen - 1; i++) {
			for (int j = 0; j < fieldLonLen; j++) {
				this->solution[i][j] = 0.5*mass->solution[i][j] * (pow(u->solution[i][j], 2) + pow(v->solution[i][j], 2));
			}
		}
	case QUADRATIC:
		for (int i = 0; i < fieldLatLen - 1; i++) {
			for (int j = 0; j < fieldLonLen; j++) {
				this->solution[i][j] = 0.5*mass->solution[i][j] * pow((pow(u->solution[i][j], 2) + pow(v->solution[i][j], 2)),1.5);
			}
		}
	}

	
};

void Energy::UpdateDtKinEAvg(void) {
	double kineticSum = 0; //Joules
	
	for (int i = 0; i < fieldLatLen-1; i++) {
		for (int j = 0; j < fieldLonLen; j++) {
			kineticSum += this->solution[i][j];
		}
	}
	//kineticSum = kineticSum*(4 / 3)*pi*(pow((consts->radius.Value() + consts->h.Value()), 3) - pow(consts->radius.Value(), 3));
	dtKinEAvg.push_back(kineticSum / (4 * pi*pow((consts->radius.Value()+consts->h.Value()),2))); //Joules per meter^2


	//Automatically update dissipation - ensures vectors of samelength
	UpdateDtDissEAvg();
};

void Energy::UpdateDtDissEAvg(void) {
	dtDissEAvg.push_back(0);

	switch (consts->fric_type) {
		//Linear dissipation
	case LINEAR:
		dtDissEAvg[dtKinEAvg.size() - 1] = 2 * dtKinEAvg[dtKinEAvg.size() - 1] * consts->alpha.Value(); //Joules per meter
	case QUADRATIC:
		dtDissEAvg[dtKinEAvg.size() - 1] = 2 * dtKinEAvg[dtKinEAvg.size() - 1] * consts->h.Value() * consts->alpha.Value();
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
		orbitDissEAvg[orbitKinEAvg.size() - 1] = 2*orbitKinEAvg[orbitKinEAvg.size() - 1]*consts->alpha.Value();
	}

};

void Energy::IsConverged(void) {
	residual.push_back(abs(orbitDissEAvg[orbitKinEAvg.size() - 1] - orbitDissEAvg[orbitKinEAvg.size() - 2]));

	if (residual.size() > 5) {
		//check latest value for convergence 
		if (residual[residual.size() - 1] < 1e-6) {
			converged = true;

			//check previous two values for convergence also:
			//Convergence will be reset if convergence is not consistent over three orbits.
			for (int i = residual.size() - 2; i > residual.size() - 5; i--) {
				if (residual[i] > 1e-6) {
					converged = false; //reset if previous two values not converged
					break;
				}
			}
		}
	}
	

	//else if (residual[residual.size() - 1] > 1) {
	//	consts->outstring << std::endl << "Residual is now greater than 1. Your model appears to have blown up. Sorry Chum." << std::endl;
	//	consts->Output.Write(ERR_MESSAGE, &consts->outstring);
	//	consts->Output.TerminateODIS();
	//}

	printf("\t Resdiual: %1.4e \n", abs(orbitDissEAvg[orbitKinEAvg.size() - 1] - orbitDissEAvg[orbitKinEAvg.size() - 2]));

	if (converged) std::cout << "Convergence criteria met." << std::endl;
};


