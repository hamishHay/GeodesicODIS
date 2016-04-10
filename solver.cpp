#ifdef _WIN32
#include <Windows.h>
#define mkdir CreateDirectory

#elif _WIN64
#include <Windows.h>

#elif __linux__
#include <unistd.h>

#else
#error "OS not supported!"
#endif

#include "solver.h"
#include "globals.h"
#include "field.h"
#include "mesh.h"
#include "mathRoutines.h"
#include "energy.h"
#include "outFiles.h"
#include "vector"

#include <math.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cstdio>
#include <sys/stat.h>
#include <errno.h>

Solver::Solver(int type, int dump, Globals * Consts, Mesh * Grid, Field * UGradLon, Field * UGradLat, Field * VelU, Field * VelV, Field * Eta, Energy * EnergyField) {
	solverType = type;
	dumpTime = dump;
	Out = &Consts->Output;

	//Assign member pointers to passed in pointers from main.cpp
	consts = Consts;
	grid = Grid;
	// dUlon = UGradLon;
	// dUlat = UGradLat;
	u = VelU;
	v = VelV;
	eta = Eta;
	energy = EnergyField;

	etaOld = new Field(grid,0,0);
	etaOldArray = etaOld->solution;

	etaNew = new Field(grid,0,0);
	etaNewArray = etaNew->solution;

	uOld = new Field(grid,0,1);
	uOldArray = uOld->solution;

	uNew = new Field(grid,0,1);
	uNewArray = uNew->solution;

	vOld = new Field(grid,1,0);
	vOldArray = vOld->solution;

	vNew = new Field(grid,1,0);
	vNewArray = vNew->solution;

	vDissTerm = new Field(grid, 1, 0);
	vDissArray = vDissTerm->solution;

	uDissTerm = new Field(grid, 0, 1);
	uDissArray = uDissTerm->solution;

	vNorthEastAvg = new Field(grid, 1, 0);
	vNEAvgArray = vNorthEastAvg->solution;

	uSouthWestAvg = new Field(grid, 0, 1);
	uSWAvgArray = uSouthWestAvg->solution;

	dUlat = new Field(grid,1,0);
	dUlatArray = dUlat->solution;

	dUlon = new Field(grid,0,1);
	dUlonArray = dUlon->solution;

	uLatLen = u->fieldLatLen;
	uLonLen = u->fieldLonLen;
	udLon = u->dLon;
	udLat = u->dLat;

	vLatLen = v->fieldLatLen;
	vLonLen = v->fieldLonLen;
	vdLon = v->dLon;
	vdLat = v->dLat;

	etaLatLen = eta->fieldLatLen;
	etaLonLen = eta->fieldLonLen;
	etadLon = eta->dLon;
	etadLat = eta->dLat;

	dt = consts->timeStep.Value();

	cosMinusB = new double[dUlat->fieldLonLen];
	cosPlusB = new double[dUlat->fieldLonLen];
	sinMinusB = new double[dUlon->fieldLonLen];
	sinPlusB = new double[dUlon->fieldLonLen];

	if (consts->potential.Value() == "ECC_RAD") tide = ECC_RAD;
	else if (consts->potential.Value() == "ECC_LIB") tide = ECC_LIB;
	else if (consts->potential.Value() == "ECC") tide = ECC;
	else if (consts->potential.Value() == "OBLIQ") tide = OBLIQ;
	else if (consts->potential.Value() == "FULL") tide = FULL;
	else {
		outstring << "No potential forcing found." << std::endl;
		Out->Write(ERR_MESSAGE, &outstring);
		Out->TerminateODIS();
	}
};

int Solver::InitialConditions(void) {
	bool action = consts->init.Value();

	outstring << "Use initial conditions: ";

	if (action) {
		outstring << "Yes." << std::endl;
		Out->Write(OUT_MESSAGE, &outstring);
		ReadInitialConditions();
		return 1;
	}

	outstring << "No." << std::endl;
	Out->Write(OUT_MESSAGE, &outstring);
	return 0;
};

void Solver::Solve() {
	outstring << "Entering solver: ";
	switch (solverType) {
	case 0:
		outstring << "Explicit\n";
		Explicit();
		break;
	default:
		outstring << "No solver type selected.\n Terminating execution.";
	}
};

void Solver::UpdatePotential() {
	switch (tide) {
	case ECC_RAD:
		UpdateEccRadPotential();
		break;

	case ECC_LIB:
		UpdateEccLibPotential();
		break;

	case ECC:
		UpdateEccPotential();
		break;

	case OBLIQ:
		UpdateObliqPotential();
		break;

	case FULL:
		UpdateFullPotential();
		break;
	}
}

inline void Solver::UpdateEccLibPotential(void) {
	double sin2Lat = 0;
	double cosLat = 0;
	double B = consts->angVel.Value()*simulationTime;
	double A = 0.25*pow(consts->angVel.Value(), 2)*pow(consts->radius.Value(), 2)*consts->e.Value();
	double Alat, Alon = 0;

	for (int j = 0; j < dUlat->fieldLonLen; j++) {
		cosMinusB[j] = cos(2 * dUlat->lon[j] - B);
		cosPlusB[j] = cos(2 * dUlat->lon[j] + B);
	}

	for (int j = 0; j < dUlon->fieldLonLen; j++) {
		sinMinusB[j] = sin(2 * dUlon->lon[j] - B);
		sinPlusB[j] = sin(2 * dUlon->lon[j] + B);
	}

	Alat = A * -1.5;
	for (int i = 0; i < dUlat->fieldLatLen; i++) {
		//lat = 2*dUlat->lat[i];
		sin2Lat = dUlat->sin2Lat[i];
		for (int j = 0; j < dUlat->fieldLonLen; j++) {
			dUlatArray[i][j] = Alat*(sin2Lat * (7 * cosMinusB[j] - cosPlusB[j])); //P22
		}
	}

	Alon = A * 3;
	for (int i = 0; i < dUlon->fieldLatLen; i++) {
		cosLat = dUlon->cosLat[i];
		for (int j = 0; j < dUlon->fieldLonLen; j++) {
			dUlonArray[i][j] = Alon * cosLat*cosLat*(-7 * sinMinusB[j] + sinPlusB[j]); //P22
		}
	}
}

inline void Solver::UpdateEccPotential(void) {
	double cosLat = 0;
	double sin2Lat = 0;
	double cosB = 0;
	double B = consts->angVel.Value()*simulationTime;
	double A = 0.25*pow(consts->angVel.Value(), 2)*pow(consts->radius.Value(), 2)*consts->e.Value();

	for (int j = 0; j < dUlat->fieldLonLen; j++) {
		cosMinusB[j] = cos(2 * dUlat->lon[j] - B);
		cosPlusB[j] = cos(2 * dUlat->lon[j] + B);
	}

	cosB = cos(B);
	for (int i = 0; i < dUlat->fieldLatLen; i++) {
		sin2Lat = dUlat->sin2Lat[i];
		for (int j = 0; j < dUlat->fieldLonLen; j++) {
			//lon = dUlat->lon[j];
			dUlatArray[i][j] = A*((-1.5*sin2Lat * (7 * cosMinusB[j] - cosPlusB[j])) + (-9.*sin2Lat*cosB)); //P22 + P20
		}
	}

	for (int j = 0; j < dUlon->fieldLonLen; j++) {
		sinMinusB[j] = sin(2 * dUlon->lon[j] - B);
		sinPlusB[j] = sin(2 * dUlon->lon[j] + B);
	}

	for (int i = 0; i < dUlon->fieldLatLen; i++) {
		cosLat = dUlon->cosLat[i];
		for (int j = 0; j < dUlon->fieldLonLen; j++) {
			//lon = dUlon->lon[j];
			dUlonArray[i][j] = A * 3 * cosLat*cosLat*(-7 * sinMinusB[j] + sinPlusB[j]); //P22
		}
	}
}

inline void Solver::UpdateEccRadPotential(void) {
	double cosB = 0;
	double sin2Lat = 0;
	double B = consts->angVel.Value()*simulationTime;
	double A = 0.25*pow(consts->angVel.Value(), 2)*pow(consts->radius.Value(), 2)*consts->e.Value();
	double value = 0;

	cosB = cos(B);
	for (int i = 0; i < dUlat->fieldLatLen; i++) {
		sin2Lat = dUlat->sin2Lat[i]; //sin(2*lat)
		value = A*(-9.*sin2Lat*cosB);
		for (int j = 0; j < dUlat->fieldLonLen; j++) {
			dUlatArray[i][j] = value; // P20
		}
	}
}

inline void Solver::UpdateObliqPotential(void) {
	double cos2Lat = 0;
	double sin2Lat = 0;
	double value = 0;
	double B = consts->angVel.Value()*simulationTime;
	double A = 0.25*pow(consts->angVel.Value(), 2)*pow(consts->radius.Value(), 2)*consts->theta.Value();

	for (int j = 0; j < dUlat->fieldLonLen; j++) {
		cosMinusB[j] = cos(dUlat->lon[j] - B);
		cosPlusB[j] = cos(dUlat->lon[j] + B);
	}

	for (int i = 0; i < dUlat->fieldLatLen; i++) {
		//lat = dUlat->lat[i];
		cos2Lat = dUlat->cos2Lat[i];
		value = -6 * A*cos2Lat;
		for (int j = 0; j < dUlat->fieldLonLen; j++) {
			//lon = dUlat->lon[j];
			dUlatArray[i][j] = value*(cosMinusB[j] + cosPlusB[j]); //P21
		}
	}

	for (int j = 0; j < dUlon->fieldLonLen; j++) {
		sinMinusB[j] = sin(dUlon->lon[j] - B);
		sinPlusB[j] = sin(dUlon->lon[j] + B);
	}

	for (int i = 0; i < dUlon->fieldLatLen; i++) {
		//lat = dUlon->lat[i]  ;
		sin2Lat = dUlon->sin2Lat[i];
		value = 3 * A * sin2Lat;
		for (int j = 0; j < dUlon->fieldLonLen; j++) {
			//lon = dUlon->lon[j];
			dUlonArray[i][j] = value*(sinMinusB[j] + sinPlusB[j]); //P21
		}
	}
}

inline void Solver::UpdateFullPotential(void) {
	double lon = 0;
	double B = consts->angVel.Value()*simulationTime;
	double A = 0.25 * pow(consts->angVel.Value(),2)*pow(consts->radius.Value(),2);
	double sin2Lat = 0;
	double cos2Lat = 0;
  double cosLat = 0;

	for (int i = 0; i < dUlat->fieldLatLen; i++) {
		sin2Lat = dUlat->sin2Lat[i];
		cos2Lat = dUlat->cos2Lat[i];
		for (int j = 0; j < dUlat->fieldLonLen; j++) {
			lon = dUlat->lon[j];
			dUlatArray[i][j] = consts->theta.Value()*-6 * cos2Lat*(cos(lon - B) + cos(lon + B));
			dUlatArray[i][j] += consts->e.Value()*(-9.*sin2Lat*cos(B));
			dUlatArray[i][j] += consts->e.Value()*(-1.5*sin2Lat * (7 * cos(2 * lon - B) - cos(2 * lon + B)));
			dUlatArray[i][j] *= A;
		}
	}

	for (int i = 0; i < dUlon->fieldLatLen; i++) {
		sin2Lat = dUlon->sin2Lat[i];
		cosLat = dUlon->cosLat[i];
		for (int j = 0; j < dUlon->fieldLonLen; j++) {
			lon = dUlon->lon[j];
			dUlonArray[i][j] = consts->theta.Value() * 3 * sin2Lat*(sin(lon - B) + sin(lon + B));
			dUlonArray[i][j] += consts->e.Value() * 3 * cosLat*cosLat*(-7 * sin(2 * lon - B) + sin(2 * lon + B));
			dUlonArray[i][j] *= A;
		}
	}
}

inline void Solver::InterpPole(Field * field) {
	for (int j = 0; j < field->fieldLonLen; j++) {
		// field->solution[0][j] = linearInterp1(field, 0, j);
		// field->solution[field->fieldLatLen - 1][j] = linearInterp1(field, field->fieldLatLen - 1, j);
		field->solution[0][j] = lagrangeInterp3(field, 0, j);
		field->solution[field->fieldLatLen - 1][j] = lagrangeInterp3(field, field->fieldLatLen - 1, j);
	}

};

inline void Solver::UpdateEastVel(){
	double coriolis = 0;
	double tidalForce = 0;

	double dSurfLon = 0;
	double surfHeight = 0;
	double eastEta = 0;
	double westEta = 0;

	double alpha = consts->alpha.Value();
	double angVel = consts->angVel.Value();

	for (int i = 1; i < uLatLen - 1; i++) {
		for (int j = 0; j < uLonLen; j++) {
			if (j < vLonLen - 1) {
				vNEAvgArray[i][j] = 0.25*(vOldArray[i][j] + vOldArray[i - 1][j] + vOldArray[i][j + 1] + vOldArray[i - 1][j + 1]);
			}
			else {
				vNEAvgArray[i][j] = 0.25*(vOldArray[i][j] + vOldArray[i - 1][j] + vOldArray[i][0] + vOldArray[i - 1][0]);
			}
		}
	}


	switch (consts->fric_type) {
	case LINEAR:
		for (int i = 1; i < uLatLen - 1; i++) {
			for (int j = 0; j < uLonLen; j++) {
				uDissArray[i][j] = uOldArray[i][j] * alpha;
			}
		}
		break;

	case QUADRATIC:
		double alphah = alpha / consts->h.Value();

		 for (int i = 1; i < uLatLen - 1; i++) {
		 	for (int j = 0; j < uLonLen; j++) {
				uDissArray[i][j] = alphah * uOldArray[i][j] * sqrt(uOldArray[i][j]*uOldArray[i][j] + vNEAvgArray[i][j]*vNEAvgArray[i][j]);
			}
		}
		break;
	}

	double loveRadius = consts->loveReduct.Value() / consts->radius.Value();
	double gRadius = consts->g.Value() / consts->radius.Value();
	double coriolisFactor = 0;
	double tidalFactor = 0;
	double surfFactor = 0;

	double * uCosLat = uOld->cosLat;
	double * uSinLat = uOld->sinLat;


	for (int i = 1; i < uLatLen - 1; i++) {
		coriolisFactor =  2. * angVel*uSinLat[i];
		tidalFactor = loveRadius/uCosLat[i];
		surfFactor = gRadius/uCosLat[i];
		for (int j = 0; j < uLonLen; j++) {
			if (j != uLonLen - 1) eastEta = etaOldArray[i][j+1];
			else eastEta = etaOldArray[i][0];

			westEta = etaOldArray[i][j];

			dSurfLon = (eastEta - westEta) / (etadLon);

			surfHeight = surfFactor*dSurfLon;

			coriolis = coriolisFactor*vNEAvgArray[i][j];
			tidalForce = tidalFactor * dUlonArray[i][j];

			uNewArray[i][j] = (coriolis - surfHeight + tidalForce - uDissArray[i][j])*dt + uOldArray[i][j];


		}
	}

	for (int j = 0; j < uLonLen; j++) {
		uNewArray[0][j] = linearInterp1Array(u,uNewArray, 0, j);
		uNewArray[uLatLen - 1][j] = linearInterp1Array(u,uNewArray, uLatLen - 1, j);
    // uNewArray[0][j] = lagrangeInterp3ArrayCenter(u,uNewArray, 0, j);
    // uNewArray[uLatLen - 1][j] = lagrangeInterp3ArrayCenter(u,uNewArray, uLatLen - 1, j);
		// uNewArray[0][j] = uNewArray[1][j];
	  // uNewArray[uLatLen - 1][j] = uNewArray[uLatLen - 2][j];
	}
	double npoleSum = 0;
	double spoleSum = 0;
	for (int j = 0; j < uLonLen; j++) {
		npoleSum += uNewArray[0][j];
		spoleSum += uNewArray[uLatLen - 1][j];
	}
	npoleSum = npoleSum / uLonLen;
	spoleSum = spoleSum / uLonLen;

	for (int j = 0; j < uLonLen; j++) {
		uNewArray[0][j] = npoleSum;
		uNewArray[uLatLen - 1][j] = spoleSum;
	}
}

inline void Solver::UpdateNorthVel(){
	double coriolis = 0;
	double tidalForce = 0;
	double dSurfLat = 0;
	double surfHeight = 0;

	double northEta = 0;
	double southEta = 0;

	double alpha = consts->alpha.Value();
	double angVel = consts->angVel.Value();

	for (int i = 0; i < vLatLen; i++) {
		for (int j = 0; j < vLonLen; j++) {
			if (j > 0) {
				uSWAvgArray[i][j] = 0.25*(uOldArray[i][j] + uOldArray[i + 1][j] + uOldArray[i][j - 1] + uOldArray[i + 1][j - 1]);
			}
			else {
				uSWAvgArray[i][j] = 0.25*(uOldArray[i][j] + uOldArray[i + 1][j] + uOldArray[i][uLonLen - 1] + uOldArray[i + 1][uLonLen - 1]);
			}
		}

	}

	switch (consts->fric_type) {
	case LINEAR:
		for (int i = 0; i < vLatLen; i++) {
			for (int j = 0; j < vLonLen; j++) {
				vDissArray[i][j] = vOldArray[i][j] * alpha;
			}
		}
		break;

	case QUADRATIC:
		double alphah = alpha / consts->h.Value();
		for (int i = 0; i < vLatLen; i++) {
			for (int j = 0; j < vLonLen; j++) {
				vDissArray[i][j] = alphah * vOldArray[i][j] * sqrt(vOldArray[i][j] * vOldArray[i][j] + uSWAvgArray[i][j]*uSWAvgArray[i][j]);
			}
		}
		break;
	}

	double loveRadius = consts->loveReduct.Value() / consts->radius.Value();
	double gRadius = consts->g.Value() / consts->radius.Value();
	double coriolisFactor = 0;

	for (int i = 0; i < vLatLen; i++) {
		coriolisFactor = 2. * angVel * v->sinLat[i];
		for (int j = 0; j < vLonLen; j++) {
			northEta = etaOldArray[i][j];
			southEta = etaOldArray[i+1][j];

			dSurfLat = (northEta - southEta) / etadLat;


			surfHeight = gRadius*dSurfLat;

			coriolis =  coriolisFactor*uSWAvgArray[i][j];

			tidalForce = loveRadius * dUlatArray[i][j];

			vNewArray[i][j] = (-coriolis - surfHeight + tidalForce - vDissArray[i][j])*dt + vOldArray[i][j];

		}
	}

}

inline void Solver::UpdateSurfaceHeight(){
	double vGrad = 0;
	double uGrad = 0;
	double northv = 0;
	double southv = 0;
	double eastu = 0;
	double westu = 0;

	double cosLat;
	double vdLat = v->dLat;
	double vdLon = u->dLon;


	double h = consts->h.Value();
	double hRadius = h / consts->radius.Value();

	for (int i = 1; i < etaLatLen-1; i++) {
		cosLat = eta->cosLat[i];
		for (int j = 0; j < etaLonLen; j++) {
			northv = vNewArray[i - 1][j] * v->cosLat[i - 1];
			southv = vNewArray[i][j] * v->cosLat[i];

			vGrad = (northv - southv) / vdLat;


			if (j != 0) {
				eastu = uNewArray[i][j];
				westu = uNewArray[i][j - 1];
			}
			else {
				eastu = uNewArray[i][j];
				westu = uNewArray[i][uLonLen - 1];
			}

			uGrad = (eastu - westu) / vdLon;

			etaNewArray[i][j] = hRadius/cosLat*(-vGrad - uGrad)*dt + etaOldArray[i][j];
		}
	}

	for (int j = 0; j < etaLonLen; j++) {
		// etaNewArray[0][j] = linearInterp1Array(eta,etaNewArray, 0, j);
		// etaNewArray[etaLatLen - 1][j] = linearInterp1Array(eta,etaNewArray, etaLatLen - 1, j);
    etaNewArray[0][j] = lagrangeInterp3ArrayCenter(eta,etaNewArray, 0, j);
    etaNewArray[etaLatLen - 1][j] = lagrangeInterp3ArrayCenter(eta,etaNewArray, etaLatLen - 1, j);
	}

	//Average eta out at poles
	double npoleSum = 0;
	double spoleSum = 0;
	for (int j = 0; j < etaLonLen; j++) {
		npoleSum += etaNewArray[0][j];
		spoleSum += etaNewArray[etaLatLen - 1][j];
	}
	npoleSum = npoleSum / etaLonLen;
	spoleSum = spoleSum / etaLonLen;

	for (int j = 0; j < etaLonLen; j++) {
		etaNewArray[0][j] = npoleSum;
		etaNewArray[etaLatLen - 1][j] = spoleSum;
	}
}


void Solver::Explicit() {
	InitialConditions();

	//Check for stability
	outstring << "Entering time loop:\n\n";
	outstring << "End time: \t" << consts->endTime.Value() / 86400.0 << " days\n";
	outstring << "Time step: \t" << dt << " seconds\n\n";
	Out->Write(OUT_MESSAGE, &outstring);

	// int output = 0;
	int outCount = 0;
	// int mun = (int) consts->period.Value()/consts->timeStep.Value()*consts->outputTime.Value();
	// int outPos[5] = {0, mun*0, mun*1, mun*2, mun*3};

	double timeStepCount = 0;
	int inc = (int) (consts->period.Value()/dt);

	DumpSolutions(-1);

	//Update cell energies and globally averaged energy
	energy->UpdateKinE(uNewArray,vNewArray);
	// energy->UpdateDtKinEAvg();

	while (simulationTime <= consts->endTime.Value() && !energy->converged) {


		timeStepCount+=dt;
		outputCount+=dt;

		UpdatePotential();

		simulationTime += dt;

		//Solve for v
		UpdateNorthVel();

		//solve for u
		UpdateEastVel();

		//Solve for eta based on new u and v
		UpdateSurfaceHeight();

		for (int i = 0; i < vLatLen; i++) {
			for (int j = 0; j < vLonLen; j++) {
				vOldArray[i][j] = vNewArray[i][j];
			}
		}

		for (int i = 0; i < uLatLen; i++) {
			for (int j = 0; j < uLonLen; j++) {
				uOldArray[i][j] = uNewArray[i][j];
			}
		}

		for (int i = 0; i < etaLatLen; i++) {
			for (int j = 0; j < etaLonLen; j++) {
				etaOldArray[i][j] = etaNewArray[i][j];
			}
		}

		energy->UpdateKinE(uNewArray,vNewArray);

		energy->UpdateDtKinEAvg();

		//Check for output
		if (timeStepCount >= consts->period.Value()) {
			orbitNumber++;
			timeStepCount -= consts->period.Value();

			energy->UpdateOrbitalKinEAvg(inc);

			// Check for convergence
			energy->IsConverged();
			if (energy->converged) convergeCount++;

			outstring << std::fixed << std::setprecision(2) << simulationTime / 86400.0 << " days: \t" << 100 * (simulationTime / consts->endTime.Value()) << "%\t" << output;
			Out->Write(OUT_MESSAGE, &outstring);

			output++;
			outCount++;
			DumpSolutions(-2);
			outCount = 1;

			energy->timePos = 0; //Reset time position after energy data output
			// DumpFields(orbitNumber);
			DumpFields(output);

		}
		else if (timeStepCount >= consts->period.Value()*consts->outputTime.Value()*outCount) {
			output++;
			outCount++;
			DumpSolutions(1);
			DumpFields(output);
		}
		iteration++;

	}

	//deallocate memory assigned to temporary fields
	delete etaNew;
	delete uNew;
	delete vNew;

	outstring << "\nSimulation complete.\n\nEnd time: \t\t" << simulationTime / 86400.0 << "\n";
	outstring << "Total iterations: \t" << iteration;
	Out->Write(OUT_MESSAGE, &outstring);

};

void Solver::DumpSolutions(int out_num) {

	if (out_num == -1) {

		//Only for windows

#if _WIN32
		mkdir(&(consts->path +  SEP + "EastVelocity" + SEP)[0], NULL);
		mkdir(&(consts->path + SEP + "NorthVelocity" + SEP)[0], NULL);
		mkdir(&(consts->path + SEP + "Displacement" + SEP)[0], NULL);
		mkdir(&(consts->path + SEP + "Grid" + SEP)[0], NULL);
		mkdir(&(consts->path + SEP + "Energy" + SEP)[0], NULL);

#elif __linux__
		mkdir(&(consts->path +  SEP + "EastVelocity" + SEP)[0], S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP);
		mkdir(&(consts->path +  SEP + "NorthVelocity" + SEP)[0], S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP);
		mkdir(&(consts->path +  SEP + "Displacement" + SEP)[0], S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP);
		mkdir(&(consts->path +  SEP + "Grid" + SEP)[0], S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP);
		mkdir(&(consts->path +  SEP + "Energy" + SEP)[0], S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP);

#endif
		//remove(&(consts->path + SEP + "diss.txt")[0]);

		std::ofstream uLat(consts->path + SEP + "Grid" + SEP + "u_lat.txt", std::ofstream::out);
		std::ofstream uLon(consts->path + SEP + "Grid" + SEP + "u_lon.txt", std::ofstream::out);
		std::ofstream vLat(consts->path + SEP + "Grid" + SEP + "v_lat.txt", std::ofstream::out);
		std::ofstream vLon(consts->path + SEP + "Grid" + SEP + "v_lon.txt", std::ofstream::out);

		for (int j = 0; j < u->ReturnFieldLonLen(); j++) uLon << u->lon[j] * 1 / radConv << '\t';
		for (int i = 0; i < u->ReturnFieldLatLen(); i++) uLat << u->lat[i] * 1 / radConv << '\t';

		for (int j = 0; j < v->ReturnFieldLonLen(); j++) vLon << v->lon[j] * 1 / radConv << '\t';
		for (int i = 0; i < v->ReturnFieldLatLen(); i++) vLat << v->lat[i] * 1 / radConv << '\t';
	}
	else if (out_num == -2) {
		FILE * dissFile;
		if (consts->kinetic.Value())
		{
			dissFile = fopen(&(consts->path + SEP + "Energy" + SEP + "kinetic_energy.txt")[0], "a+");
			fprintf(dissFile, "%.30f \n", energy->dtKinEAvg[energy->timePos-1]);
			fclose(dissFile);

			dissFile = fopen(&(consts->path + SEP + "Energy" + SEP + "kinetic_energy_orb_avg.txt")[0], "a+");
			fprintf(dissFile, "%.30f \n", energy->orbitKinEAvg);
			fclose(dissFile);
		}
		if (consts->diss.Value()) {
			dissFile = fopen(&(consts->path + SEP + "Energy" + SEP + "diss_energy.txt")[0], "a+");
			fprintf(dissFile, "%.30f \n", energy->dtDissEAvg[energy->timePos-1]);
			fclose(dissFile);

			dissFile = fopen(&(consts->path + SEP + "Energy" + SEP + "diss_energy_orb_avg.txt")[0], "a+");
			fprintf(dissFile, "%.30f \n", energy->orbitDissEAvg[0]);
			fclose(dissFile);
		}
	}
	else {
		//fopen for linux
		//fopen_s for windows
		FILE * dissFile;
		if (consts->kinetic.Value())
		{
			dissFile = fopen(&(consts->path + SEP + "Energy" + SEP + "kinetic_energy.txt")[0], "a+");
			fprintf(dissFile, "%.30f \n", energy->dtKinEAvg[energy->timePos-1]);
			fclose(dissFile);
		}
		if (consts->diss.Value())
		{
			dissFile = fopen(&(consts->path + SEP + "Energy" + SEP + "diss_energy.txt")[0], "a+");
			fprintf(dissFile, "%.30f \n", energy->dtDissEAvg[energy->timePos-1]);
			fclose(dissFile);
		}
		if (consts->work.Value())
		{
			dissFile = fopen(&(consts->path + SEP + "Energy" + SEP + "energy.txt")[0], "a+");
			// fprintf(dissFile, "%.15f \n", energy->orbitDissEAvg[1]);
			fprintf(dissFile, "%.30f \n", energy->orbitKinEAvg);
			// fprintf(dissFile, "%.15f \n", energy->dtKinEAvg[energy->timePos-1]);
			// fprintf(dissFile, "%.15f \n", energy->dtDissEAvg[energy->timePos-1]);
			// fprintf(dissFile, "%.15f \n", energy->dissipation[energy->dissipation.size()-1]);
			fclose(dissFile);
		}

		if (energy->converged) DumpFields(out_num);
		// else if (out_num % 5 == 0) DumpFields(out_num); //dump every 5 orbits

		// DumpFields(out_num); //FOR DEBUGGING
	}
};

void Solver::ReadInitialConditions(void) {
	std::ifstream eastVel(consts->path + SEP + "InitialConditions" + SEP + "u_vel.txt", std::ifstream::in);
	std::ifstream northVel(consts->path + SEP + "InitialConditions" + SEP + "v_vel.txt", std::ifstream::in);
	std::ifstream displacement(consts->path + SEP + "InitialConditions" + SEP + "eta.txt", std::ifstream::in);

	CopyInitialConditions(eastVel, uOld);
	CopyInitialConditions(northVel, vOld);
	CopyInitialConditions(displacement, etaOld);
};

void Solver::CopyInitialConditions(std::ifstream & file, Field * inField) {
	double inputValue = 0;
	std::string line;

	file.is_open();

	for (int i = 0; i < inField->ReturnFieldLatLen(); i++) {
		for (int j = 0; j < inField->ReturnFieldLonLen(); j++) {
			std::getline(file >> std::ws, line, '\t');
			std::stringstream inputString(line);
			inputString >> inputValue;

			inField->solution[i][j] = inputValue;
		}
	}

	file.close();
};

void Solver::DumpFields(int output_num) {
	std::string out = std::to_string(output_num);

	std::ofstream uFile(consts->path + SEP + "EastVelocity" + SEP + "u_vel_" + out + ".txt", std::ofstream::out);
	std::ofstream vFile(consts->path + SEP + "NorthVelocity" + SEP + "v_vel_" + out + ".txt", std::ofstream::out);
	std::ofstream etaFile(consts->path + SEP + "Displacement" + SEP + "eta_" + out + ".txt", std::ofstream::out);

	for (int i = 0; i < u->ReturnFieldLatLen(); i++) {
		for (int j = 0; j < u->ReturnFieldLonLen(); j++) {
			uFile << uNewArray[i][j] << '\t';
			// etaFile << etaNew->solution[i][j] << '\t';
			etaFile << etaNewArray[i][j] << '\t';
		}
		uFile << std::endl;
		etaFile << std::endl;
	}

	for (int i = 0; i < vLatLen; i++) {
		for (int j = 0; j < vLonLen; j++) {
			vFile << vOldArray[i][j] << '\t';
		}
		vFile << std::endl;
	}
};
