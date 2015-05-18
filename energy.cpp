#include "field.h"
#include "energy.h"
#include <math.h>

KineticEnergy::KineticEnergy(Mesh * mesh, int lat, int lon) : Field (mesh, lat, lon) {};

void KineticEnergy::UpdateKinE(Field * u, Field * v) {
	double cellMass = 0.;

	for (int i = 0; i < fieldLatLen; i++) {
		for (int j = 0; j < fieldLonLen; j++) {
			this->solution[i][j] = 0.5*cellMass*(pow(u->solution[i][j], 2) + pow(v->solution[i][j], 2));
		}
	}
};

void KineticEnergy::TimeStepKinEAvg(void) {
	double kineticSum = 0;

	//Sum for poles
	//kineticSum += this->solution[0][0];
	//kineticSum += this->solution[this->fieldLatLen-1][0];
	
	for (int i = 1; i < fieldLatLen-1; i++) {
		for (int j = 0; j < fieldLonLen; j++) {
			kineticSum += this->solution[i][j];
		}
	}

	timeStepKinAvg.push_back(kineticSum / (fieldLatLen*fieldLonLen));
};

void KineticEnergy::OrbitalKinEAvg(int inc) {	
	orbitKinAvg.push_back(0);

	int pos = orbitKinAvg.size();
	
	for (int i = timeStepKinAvg.size()-1; i = timeStepKinAvg.size() - inc; i--) {
		orbitKinAvg[pos] += timeStepKinAvg[i];
	}

	orbitKinAvg[pos] /= inc;

}