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
	dUlon = UGradLon;
	dUlat = UGradLat;
	u = VelU;
	v = VelV;
	eta = Eta;
	energy = EnergyField;

	etaNew = new Field(grid,0,0);
	uNew = new Field(grid,0,1);
	vNew = new Field(grid,1,0);

	vDissTerm = new Field(grid, 1, 0);
	uDissTerm = new Field(grid, 0, 1);

	vNorthEastAvg = new Field(grid, 1, 0);
	uSouthWestAvg = new Field(grid, 0, 1);

	*etaNew = *eta;
	*uNew = *u;
	*etaNew = *eta;

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
	double lat = 0;
	double lon = 0;
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
			dUlat->solution[i][j] = Alat*(sin2Lat * (7 * cosMinusB[j] - cosPlusB[j])); //P22
		}
	}

	Alon = A * 3;
	for (int i = 0; i < dUlon->fieldLatLen; i++) {
		cosLat = dUlon->cosLat[i];
		for (int j = 0; j < dUlon->fieldLonLen; j++) {
			dUlon->solution[i][j] = Alon * cosLat*cosLat*(-7 * sinMinusB[j] + sinPlusB[j]); //P22
		}
	}
}

inline void Solver::UpdateEccPotential(void) {
	double lat = 0;
	double lon = 0;
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
			dUlat->solution[i][j] = A*((-1.5*sin2Lat * (7 * cosMinusB[j] - cosPlusB[j])) + (-9.*sin2Lat*cosB)); //P22 + P20
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
			dUlon->solution[i][j] = A * 3 * cosLat*cosLat*(-7 * sinMinusB[j] + sinPlusB[j]); //P22
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
			dUlat->solution[i][j] = value; // P20
		}
	}
}

inline void Solver::UpdateObliqPotential(void) {
	//double lat = 0;
	//double lon = 0;
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
			dUlat->solution[i][j] = value*(cosMinusB[j] + cosPlusB[j]); //P21
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
			dUlon->solution[i][j] = value*(sinMinusB[j] + sinPlusB[j]); //P21
		}
	}
}

inline void Solver::UpdateFullPotential(void) {
	double lat = 0;
	double lon = 0;
	double B = consts->angVel.Value()*simulationTime;
	double A = 0.25 * pow(consts->angVel.Value(),2)*pow(consts->radius.Value(),2);

	for (int i = 0; i < dUlat->fieldLatLen; i++) {
		lat = dUlat->lat[i]  ;
		for (int j = 0; j < dUlat->fieldLonLen; j++) {
			lon = dUlat->lon[j];
			dUlat->solution[i][j] = consts->theta.Value()*-6 * cos(2 * lat)*(cos(lon - B) + cos(lon + B));
			dUlat->solution[i][j] += consts->e.Value()*(-9.*sin(2 * lat)*cos(B));
			dUlat->solution[i][j] += consts->e.Value()*(-1.5*sin(2 * lat) * (7 * cos(2 * lon - B) - cos(2 * lon + B)));
			dUlat->solution[i][j] *= A;
		}
	}

	for (int i = 0; i < dUlon->fieldLatLen; i++) {
		lat = dUlon->lat[i]  ;
		for (int j = 0; j < dUlon->fieldLonLen; j++) {
			lon = dUlon->lon[j]  ;
			dUlon->solution[i][j] = consts->theta.Value() * 3 * sin(2 * lat)*(sin(lon - B) + sin(lon + B));
			dUlon->solution[i][j] += consts->e.Value() * 3 * cos(lat)*cos(lat)*(-7 * sin(2 * lon - B) + sin(2 * lon + B));
			dUlon->solution[i][j] *= A;
		}
	}
}

inline void Solver::InterpPole(Field * field) {
	//double L[4];
	//double x[5];
	//int inc = 1;

	//// North pole:
	//int i = 0;
	//for (int k = 0; k < 2; k++) {
	//	for (int j = 0; j < field->fieldLonLen; j++) {

	//		x[0] = field->lat[i];
	//		x[1] = field->lat[i + inc];
	//		x[2] = field->lat[i + inc * 2];
	//		x[3] = field->lat[i + inc * 3];
	//		x[4] = field->lat[i + inc * 4];

	//		//solve for each quadratic term
	//		L[0] = (x[0] - x[2])*(x[0] - x[3])*(x[0] - x[4]) / ((x[1] - x[2])*(x[1] - x[3])*(x[1] - x[4])) * field->solution[i + inc][j]; //L1
	//		L[1] = (x[0] - x[1])*(x[0] - x[3])*(x[0] - x[4]) / ((x[2] - x[1])*(x[2] - x[3])*(x[2] - x[4])) * field->solution[i + inc * 2][j]; //L2
	//		L[2] = (x[0] - x[1])*(x[0] - x[2])*(x[0] - x[4]) / ((x[3] - x[1])*(x[3] - x[2])*(x[3] - x[4])) * field->solution[i + inc * 3][j]; //L3
	//		L[3] = (x[0] - x[1])*(x[0] - x[2])*(x[0] - x[3]) / ((x[4] - x[1])*(x[4] - x[2])*(x[4] - x[3])) * field->solution[i + inc * 4][j];

	//		field->solution[i][j] = L[0] + L[1] + L[2] + L[3];
	//	}

	//	inc = -inc;
	//	i = field->fieldLatLen - 1;

	//}
	for (int j = 0; j < field->fieldLonLen; j++) {
		field->solution[0][j] = linearInterp1(field, 0, j);
		field->solution[field->fieldLatLen - 1][j] = linearInterp1(field, field->fieldLatLen - 1, j);
	}

};

inline void Solver::UpdateEastVel(Field * U, Field * UNEW, Field * V, Field * ETA){
	double coriolis = 0;
	double tidalForce = 0;

	double dSurfLon = 0;
	double surfHeight = 0;
	double eastEta = 0;
	double westEta = 0;

	double northEastAvg = 0;

//	double cosLat = 0;
//	double sinLat = 0;
	double alpha = consts->alpha.Value();
	double angVel = consts->angVel.Value();

	for (int i = 1; i < u->fieldLatLen - 1; i++) {
		for (int j = 0; j < u->fieldLonLen; j++) {
			if (j < V->fieldLonLen - 1) {
				vNorthEastAvg->solution[i][j] = 0.25*(V->solution[i][j] + V->solution[i - 1][j] + V->solution[i][j + 1] + V->solution[i - 1][j + 1]);
			}
			else {
				vNorthEastAvg->solution[i][j] = 0.25*(V->solution[i][j] + V->solution[i - 1][j] + V->solution[i][0] + V->solution[i - 1][0]);
			}
		}
	}

	switch (consts->fric_type) {
	case LINEAR:
		for (int i = 1; i < u->fieldLatLen - 1; i++) {
			for (int j = 0; j < u->fieldLonLen; j++) {
				uDissTerm->solution[i][j] = U->solution[i][j] * alpha;
			}
		}
		break;

	case QUADRATIC:
		double alphah = alpha / consts->h.Value();

		 for (int i = 1; i < u->fieldLatLen - 1; i++) {
		 	for (int j = 0; j < u->fieldLonLen; j++) {
				northEastAvg = vNorthEastAvg->solution[i][j];
				uDissTerm->solution[i][j] = alphah * U->solution[i][j] * sqrt(U->solution[i][j]*U->solution[i][j] + northEastAvg*northEastAvg);
			}
		}
		break;
	}

	double loveRadius = consts->loveReduct.Value() / consts->radius.Value();
	double gRadius = consts->g.Value() / consts->radius.Value();
	double coriolisFactor = 0;
	double tidalFactor = 0;
	double surfFactor = 0;


	for (int i = 1; i < u->fieldLatLen - 1; i++) {
		//cosLat = u->cosLat[i];
		coriolisFactor =  2. * angVel*u->sinLat[i];
		tidalFactor = loveRadius/u->cosLat[i];
		surfFactor = gRadius/u->cosLat[i];
		for (int j = 0; j < u->fieldLonLen; j++) {
			if (j != u->fieldLonLen - 1) eastEta = ETA->solution[i][j + 1];
			else eastEta = ETA->solution[i][0];

			westEta = ETA->solution[i][j];

			dSurfLon = (eastEta - westEta) / (eta->dLon);

			surfHeight = surfFactor*dSurfLon;

			coriolis = coriolisFactor*vNorthEastAvg->solution[i][j];
			tidalForce = tidalFactor * dUlon->solution[i][j];

			UNEW->solution[i][j] = (coriolis - surfHeight + tidalForce - uDissTerm->solution[i][j])*dt + U->solution[i][j];
		}
	}


	InterpPole(UNEW);
}

inline void Solver::UpdateNorthVel(Field * V, Field * VNEW, Field * U, Field * ETA){
	double coriolis = 0;
	double tidalForce = 0;
	double dSurfLat = 0;
	double surfHeight = 0;

	double northEta = 0;
	double southEta = 0;
	double etadLat = eta->dLat;

	double southwestavg = 0;

	double sinLat;

	double alpha = consts->alpha.Value();
	double angVel = consts->angVel.Value();

	for (int i = 0; i < v->fieldLatLen; i++) {
		for (int j = 0; j < v->fieldLonLen; j++) {
			if (j > 0) {
				uSouthWestAvg->solution[i][j] = 0.25*(U->solution[i][j] + U->solution[i + 1][j] + U->solution[i][j - 1] + U->solution[i + 1][j - 1]);
			}
			else {
				uSouthWestAvg->solution[i][j] = 0.25*(U->solution[i][j] + U->solution[i + 1][j] + U->solution[i][U->fieldLonLen - 1] + U->solution[i + 1][U->fieldLonLen - 1]);
			}
		}
	}

	switch (consts->fric_type) {
	case LINEAR:
		for (int i = 0; i < v->fieldLatLen; i++) {
			for (int j = 0; j < v->fieldLonLen; j++) {
				vDissTerm->solution[i][j] = V->solution[i][j] * alpha;
			}
		}
		break;

	case QUADRATIC:
		double alphah = alpha / consts->h.Value();
		for (int i = 0; i < v->fieldLatLen; i++) {
			for (int j = 0; j < v->fieldLonLen; j++) {
				southwestavg = uSouthWestAvg->solution[i][j];
				vDissTerm->solution[i][j] = alphah * V->solution[i][j] * sqrt(V->solution[i][j] * V->solution[i][j] + southwestavg*southwestavg);
			}
		}
		break;
	}

	double loveRadius = consts->loveReduct.Value() / consts->radius.Value();
	double gRadius = consts->g.Value() / consts->radius.Value();
	double coriolisFactor = 0;

	for (int i = 0; i < v->fieldLatLen; i++) {
		coriolisFactor = 2. * angVel * v->sinLat[i];
		for (int j = 0; j < v->fieldLonLen; j++) {
			northEta = eta->solution[i][j];
			southEta = eta->solution[i + 1][j];

			dSurfLat = (northEta - southEta) / etadLat;

			surfHeight = gRadius*dSurfLat;

			coriolis =  coriolisFactor*uSouthWestAvg->solution[i][j];

			tidalForce = loveRadius * dUlat->solution[i][j];

			VNEW->solution[i][j] = (-coriolis - surfHeight + tidalForce - vDissTerm->solution[i][j])*dt + V->solution[i][j];
		}
	}
}

inline void Solver::UpdateSurfaceHeight(Field * ETA, Field * ETANEW, Field * U, Field * V){
	double vGrad = 0;
	double uGrad = 0;
	double northv = 0;
	double southv = 0;
	double eastu = 0;
	double westu = 0;

	double cosLat;
	double vdLat = v->dLat;
	double vdLon = u->dLon;

	double hRadius = consts->h.Value() / consts->radius.Value();

	for (int i = 1; i < eta->fieldLatLen-1; i++) {
		cosLat = eta->cosLat[i];
		for (int j = 0; j < eta->fieldLonLen; j++) {
			northv = V->solution[i - 1][j] * v->cosLat[i - 1];
			southv = V->solution[i][j] * v->cosLat[i];

			vGrad = (northv - southv) / vdLat;

			if (j != 0) {
				eastu = U->solution[i][j];
				westu = U->solution[i][j - 1];
			}
			else {
				westu = U->solution[i][U->fieldLonLen - 1];
				eastu = U->solution[i][j];
			}

			uGrad = (eastu - westu) / vdLon;

			ETANEW->solution[i][j] = hRadius/cosLat*(-vGrad - uGrad)*dt + ETA->solution[i][j];
		}
	}

	InterpPole(ETANEW);

	//Average eta out at poles

	double npoleSum = 0;
	double spoleSum = 0;
	for (int j = 0; j < eta->fieldLonLen; j++) {
		npoleSum += ETANEW->solution[0][j];
		spoleSum += ETANEW->solution[eta->fieldLatLen - 1][j];
	}
	npoleSum = npoleSum / eta->fieldLonLen;
	spoleSum = spoleSum / eta->fieldLonLen;

	for (int j = 0; j < eta->fieldLonLen; j++) {
		ETANEW->solution[0][j] = npoleSum;
		ETANEW->solution[eta->fieldLatLen - 1][j] = spoleSum;
	}
}


void Solver::Explicit() {
	InitialConditions();

	//Check for stability
	outstring << "Entering time loop:\n\n";
	outstring << "End time: \t" << consts->endTime.Value() / 86400.0 << " days\n";
	outstring << "Time step: \t" << dt << " seconds\n\n";
	Out->Write(OUT_MESSAGE, &outstring);

	int output = 0;

	double timeStepCount = 0;
	int inc = (int) (consts->period.Value()/dt);
	DumpSolutions(-1);

	//Update cell energies and globally averaged energy
	energy->UpdateKinE();

	while (simulationTime <= consts->endTime.Value() && convergeCount <= convergeMax) {


		timeStepCount+=dt;

		UpdatePotential();

		simulationTime += dt;

		//Solve for v
		UpdateNorthVel(v, vNew, u, eta);

		//solve for u
		UpdateEastVel(u, uNew, v, eta);

		//Solve for eta based on new u and v
		UpdateSurfaceHeight(eta, etaNew, uNew, vNew);

		//Overwite previous timestep solutions at end of iteration.
		*eta = *etaNew;
		*u = *uNew;
		*v = *vNew;

		energy->UpdateKinE();

		energy->UpdateDtKinEAvg();

		//Check for output

		if (timeStepCount >= consts->period.Value()) {
		//if (true) {
			orbitNumber++;
			timeStepCount = timeStepCount - consts->period.Value();

			energy->UpdateOrbitalKinEAvg(inc);

			//Check for convergence
			if (orbitNumber>1) energy->IsConverged();
			if (energy->converged) convergeCount++;

			outstring << std::fixed << std::setprecision(2) << simulationTime / 86400.0 << " days: \t" << 100 * (simulationTime / consts->endTime.Value()) << "%\t" << output;
			Out->Write(OUT_MESSAGE, &outstring);

			DumpSolutions(output);
			output++;
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

#elif __linux__
		mkdir(&(consts->path +  SEP + "EastVelocity" + SEP)[0], S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP);
		mkdir(&(consts->path +  SEP + "NorthVelocity" + SEP)[0], S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP);
		mkdir(&(consts->path +  SEP + "Displacement" + SEP)[0], S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP);
		mkdir(&(consts->path +  SEP + "Grid" + SEP)[0], S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP);

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

	else {
		//fopen for linux
		//fopen_s for windows
		FILE * dissFile = fopen(&(consts->path + SEP + "diss.txt")[0], "a+");
		fprintf(dissFile, "%.15f \n", energy->orbitDissEAvg[1]);
		fclose(dissFile);

		if (energy->converged) DumpFields(out_num);
		else if (out_num % 20 == 0) DumpFields(out_num); //dump every 5 orbits
	}
};

void Solver::ReadInitialConditions(void) {
	std::ifstream eastVel(consts->path + SEP + "InitialConditions" + SEP + "u_vel.txt", std::ifstream::in);
	std::ifstream northVel(consts->path + SEP + "InitialConditions" + SEP + "v_vel.txt", std::ifstream::in);
	std::ifstream displacement(consts->path + SEP + "InitialConditions" + SEP + "eta.txt", std::ifstream::in);

	CopyInitialConditions(eastVel, u);
	CopyInitialConditions(northVel, v);
	CopyInitialConditions(displacement, eta);
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
			uFile << uNew->solution[i][j] << '\t';
			etaFile << etaNew->solution[i][j] << '\t';
		}
		uFile << std::endl;
		etaFile << std::endl;
	}

	for (int i = 0; i < v->fieldLatLen; i++) {
		for (int j = 0; j < v->fieldLonLen; j++) {
			vFile << vNew->solution[i][j] << '\t';
		}
		vFile << std::endl;
	}
};
