#include "field.h"
#include "energy.h"
#include "mass.h"
#include "globals.h"
#include "outFiles.h"
#include "mass.h"

#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>

Energy::Energy(Mesh * mesh, int lat, int lon, Globals * Consts, Field * UVel, Field * VVel, Mass * MassField) : Field (mesh, lat, lon)
{
	consts = Consts;
	u = UVel;
	v = VVel;
	mass = MassField;

	timePos = 0;
	totalSize = (int) (consts->period.Value()/consts->timeStep.Value());

	dtKinEAvg = new double[totalSize+1];
	dtDissEAvg = new double[totalSize+1];
};

void Energy::UpdateKinE(void) {
	switch (consts->fric_type) {
		//Linear dissipation
	case LINEAR:
		for (int i = 0; i < fieldLatLen - 1; i++) {
			for (int j = 0; j < fieldLonLen; j++) {
				this->solution[i][j] = 0.5*mass->solution[i][j] * (pow(u->solution[i][j], 2) + pow(v->solution[i][j], 2));
			}
		}
		break;

	case QUADRATIC:
		for (int i = 0; i < fieldLatLen - 1; i++) {
			for (int j = 0; j < fieldLonLen; j++) {
				this->solution[i][j] = 0.5*mass->solution[i][j] * pow((pow(u->solution[i][j], 2) + pow(v->solution[i][j], 2)), 1.5);
			}
		}
		break;
	}
};

void Energy::UpdateDtKinEAvg(void) {
	double kineticSum = 0; //Joules

	for (int i = 0; i < fieldLatLen-1; i++) {
		for (int j = 0; j < fieldLonLen; j++) {
			kineticSum += this->solution[i][j];
		}
	}
	dtKinEAvg[timePos] = kineticSum / (4 * pi*pow(consts->radius.Value(),2)); //Joules per meter^2

	//Automatically update dissipation - ensures vectors of samelength
	UpdateDtDissEAvg();
};

void Energy::UpdateDtDissEAvg(void) {
	switch (consts->fric_type) {
		//Linear dissipation
	case LINEAR:
		dtDissEAvg[timePos] = 2 * dtKinEAvg[timePos] * consts->alpha.Value(); //Joules per meter^2
		break;
	case QUADRATIC:
		dtDissEAvg[timePos] = 2 * dtKinEAvg[timePos] * consts->alpha.Value()/consts->h.Value();
		break;
	default:
		consts->Output.TerminateODIS();
	}

	timePos += 1;
};


void Energy::UpdateOrbitalKinEAvg(int inc) {
	orbitKinEAvg = 0;

	//int pos = orbitKinEAvg.size()-1;

	for (unsigned int i = 0; i < timePos; i++) {
		orbitKinEAvg += dtKinEAvg[i];
	}

	orbitKinEAvg /= inc;

	//Automatically update dissipation
	UpdateOrbitalDissEAvg();

};

void Energy::UpdateOrbitalDissEAvg(void) {
	orbitDissEAvg[1] = orbitDissEAvg[0];
	orbitDissEAvg[0] = 0;

	switch (consts->fric_type) {
		//Linear dissipation
	case LINEAR:
		orbitDissEAvg[0] = 2 * orbitKinEAvg * consts->alpha.Value(); //Joules per meter
		break;
	case QUADRATIC:
		orbitDissEAvg[0] = 2 * orbitKinEAvg * consts->alpha.Value() / consts->h.Value();
		break;
	}

	//Reset time position after updating orbital values
	timePos = 0;

	dissipation.push_back(orbitDissEAvg[0]);

};

void Energy::IsConverged(void) {
	if (dissipation.size() > 10 && !converged) {
		for (int i=derivative.size(); i < dissipation.size() - 1; i++) {
			derivative.push_back(dissipation[i+1]-dissipation[i]);
		}

		for (int i=derivative.size()-2; i < derivative.size()-1; i++) {
			// Check for minima
			if (derivative[i] < 0 && derivative[i+1] > 0) {
				minima.push_back(i);
				std::cout<<"Minima found at orbit "<<i<<std::endl;
			}
			else if (derivative[i] > 0 && derivative[i+1] < 0) {
				maxima.push_back(i);
				std::cout<<"Maxima found at orbit "<<i<<std::endl;
			}
		}

		if (maxima.size() == 2 || minima.size() == 2) converged = true;
	}

	residual.push_back(fabs(orbitDissEAvg[1] - orbitDissEAvg[0]));
	if (residual.size() > 8 && !converged) {
		//check latest value for convergence
		if (residual[residual.size() - 1] < consts->converge.Value()) {
			converged = true;

			//check previous two values for convergence also:
			//Convergence will be reset if convergence is not consistent over three orbits.
			for (unsigned int i = residual.size() - 2; i > residual.size() - 8; i--) {
				if (residual[i] > consts->converge.Value()) {
					converged = false; //reset if previous two values not converged
					break;
				}
			}
		}
	}

	if (residual[residual.size() - 1] > 1e5) {
		consts->outstring << std::endl << "Residual is now greater than 100,000. Your model appears to have blown up. Sorry Chum." << std::endl;
		consts->Output.Write(ERR_MESSAGE, &consts->outstring);
		consts->Output.TerminateODIS();
	}

	printf("\t Resdiual: %1.4e \n", residual[residual.size() - 1]);
	if (converged) std::cout << "Convergence criteria met." << std::endl;

	count += 1;
};
