#include "solver.h"
#include "globals.h"
#include "field.h"
#include "mesh.h"
#include "mathRoutines.h"
#include "mass.h"
#include "energy.h"
#include "outFiles.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <algorithm>
#include <cstdio>
#include <Windows.h>

Solver::Solver(int type, int dump, Globals * Consts, Mesh * Grid, Field * UGradLon, Field * UGradLat, Field * VelU, Field * VelV, Field * Eta, Mass * MassField, Energy * EnergyField) {
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
	mass = MassField;
	energy = EnergyField;

	etaNew = new Field(grid,0,0);
	uNew = new Field(grid,0,1);
	vNew = new Field(grid,1,0);

	*etaNew = *eta;
	*uNew = *u;
	*etaNew = *eta;

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

void Solver::UpdateEccLibPotential(void) {
	double lat = 0;
	double lon = 0;
	double B = consts->angVel.Value()*simulationTime;
	double constant = 0.25*pow(consts->angVel.Value(), 2)*pow(consts->radius.Value(), 2)*consts->e.Value();

	for (int i = 0; i < dUlat->fieldLatLen; i++) {
		lat = dUlat->lat[i] * radConv;
		for (int j = 0; j < dUlat->fieldLonLen; j++) {
			lon = dUlat->lon[j] * radConv;
			dUlat->solution[i][j] = constant*(- 1.5*sin(2*lat) * (7 * cos(2*lon - B) - cos(2*lon + B))); //P22 
		}
	}

	for (int i = 0; i < dUlon->fieldLatLen; i++) {
		lat = dUlon->lat[i] * radConv;
		for (int j = 0; j < dUlon->fieldLonLen; j++) {
			lon = dUlon->lon[j] * radConv;
			dUlon->solution[i][j] = constant *3* pow(cos(lat), 2)*(-7 * sin(2*lon - B) + sin(2*lon + B)); //P22
		}
	}
}

void Solver::UpdateEccPotential(void) {
	double lat = 0;
	double lon = 0;
	double B = consts->angVel.Value()*simulationTime;
	double constant = 0.25*pow(consts->angVel.Value(), 2)*pow(consts->radius.Value(), 2)*consts->e.Value();

	for (int i = 0; i < dUlat->fieldLatLen; i++) {
		lat = dUlat->lat[i] * radConv;
		for (int j = 0; j < dUlat->fieldLonLen; j++) {
			lon = dUlat->lon[j] * radConv;
			dUlat->solution[i][j] = constant*(-1.5*sin(2 * lat) * (7 * cos(2 * lon - B) - cos(2 * lon + B))) + constant*(-9.*sin(2 * lat)*cos(B)); //P22 + P20
		}
	}

	for (int i = 0; i < dUlon->fieldLatLen; i++) {
		lat = dUlon->lat[i] * radConv;
		for (int j = 0; j < dUlon->fieldLonLen; j++) {
			lon = dUlon->lon[j] * radConv;
			dUlon->solution[i][j] = constant * 3 * pow(cos(lat), 2)*(-7 * sin(2 * lon - B) + sin(2 * lon + B)); //P22
		}
	}
}

void Solver::UpdateEccRadPotential(void) {
	double lat = 0;
	double lon = 0;
	double B = consts->angVel.Value()*simulationTime;
	double constant = 0.25*pow(consts->angVel.Value(), 2)*pow(consts->radius.Value(), 2)*consts->e.Value();//consts->theta.Value();//

	for (int i = 0; i < dUlat->fieldLatLen; i++) {
		lat = dUlat->lat[i] * radConv;
		for (int j = 0; j < dUlat->fieldLonLen; j++) {
			lon = dUlat->lon[j] * radConv;
			dUlat->solution[i][j] = constant*(-9.*sin(2*lat)*cos(B)); // P20
		}
	}
}

void Solver::UpdateObliqPotential(void) {
	double lat = 0;
	double lon = 0;
	double B = consts->angVel.Value()*simulationTime;
	double constant = 0.25*pow(consts->angVel.Value(), 2)*pow(consts->radius.Value(), 2)*consts->theta.Value();

	for (int i = 0; i < dUlat->fieldLatLen; i++) {
		lat = dUlat->lat[i] * radConv;
		for (int j = 0; j < dUlat->fieldLonLen; j++) {
			lon = dUlat->lon[j] * radConv;
			dUlat->solution[i][j] = -6 * constant*cos(2 * lat)*(cos(lon - B) + cos(lon + B)); //P21
		}
	}

	for (int i = 0; i < dUlon->fieldLatLen; i++) {
		lat = dUlon->lat[i] * radConv;
		for (int j = 0; j < dUlon->fieldLonLen; j++) {
			lon = dUlon->lon[j] * radConv;
			dUlon->solution[i][j] = 3 * constant * sin(2 * lat)*(sin(lon - B) + sin(lon + B)); //P21
		}
	}
}

void Solver::UpdateFullPotential(void) {
	double lat = 0;
	double lon = 0;
	double B = consts->angVel.Value()*simulationTime;

	for (int i = 0; i < dUlat->fieldLatLen; i++) {
		lat = dUlat->lat[i] * radConv;
		for (int j = 0; j < dUlat->fieldLonLen; j++) {
			lon = dUlat->lon[j] * radConv;
			dUlat->solution[i][j] = 0.25*pow(consts->angVel.Value(), 2)*pow(consts->radius.Value(), 2)*consts->theta.Value()*-6 * cos(2 * lat)*(cos(lon - B) + cos(lon + B));
			dUlat->solution[i][j] += 0.25*pow(consts->angVel.Value(), 2)*pow(consts->radius.Value(), 2)*consts->e.Value()*(-9.*sin(2 * lat)*cos(B));
			dUlat->solution[i][j] += 0.25*pow(consts->angVel.Value(), 2)*pow(consts->radius.Value(), 2)*consts->e.Value()*(-1.5*sin(2 * lat) * (7 * cos(2 * lon - B) - cos(2 * lon + B)));
		}
	}

	for (int i = 0; i < dUlon->fieldLatLen; i++) {
		lat = dUlon->lat[i] * radConv;
		for (int j = 0; j < dUlon->fieldLonLen; j++) {
			lon = dUlon->lon[j] * radConv;
			dUlon->solution[i][j] = 0.25*pow(consts->angVel.Value(), 2)*pow(consts->radius.Value(), 2)*consts->theta.Value() * 3 * sin(2 * lat)*(sin(lon - B) + sin(lon + B));
			dUlon->solution[i][j] += 0.25*pow(consts->angVel.Value(), 2)*pow(consts->radius.Value(), 2)*consts->e.Value() * 3 * pow(cos(lat), 2)*(-7 * sin(2 * lon - B) + sin(2 * lon + B));
		}
	}
}

void Solver::UpdateEastVel(Field * UOLD, Field * UNEW, Field * U, Field * V, Field * ETA, double dt){
	double coriolis = 0;
	double tidalForce = 0;

	double dSurfLon = 0;
	double surfHeight = 0;
	double eastEta = 0;
	double westEta = 0;

	double lat, lon;

	for (int i = 1; i < u->fieldLatLen - 1; i++) {
		lat = u->lat[i] * radConv;
		for (int j = 0; j < u->fieldLonLen; j++) {
			lon = u->lon[j] * radConv;

			eastEta = ETA->EastP(i, j);
			westEta = ETA->CenterP(i, j);

			dSurfLon = (eastEta - westEta) / (eta->dLon*radConv);

			surfHeight = (consts->g.Value() / (consts->radius.Value()*cos(lat)))*dSurfLon;

			coriolis = 2. * consts->angVel.Value()*sin(lat)*V->NorthEastAvg(i, j);
			tidalForce = consts->loveReduct.Value()*(1. / (consts->radius.Value()*cos(lat)))* dUlon->solution[i][j];

			UNEW->solution[i][j] = (coriolis - surfHeight + tidalForce - U->solution[i][j] * consts->alpha.Value())*dt + UOLD->solution[i][j];
		}
	}

	for (int j = 0; j < eta->fieldLonLen; j++) {
		UNEW->solution[0][j] = linearInterp1(UNEW, 0, j);
		UNEW->solution[u->fieldLatLen - 1][j] = linearInterp1(UNEW, u->fieldLatLen - 1, j);
	}

}

void Solver::UpdateNorthVel(Field * VOLD, Field * VNEW, Field * U, Field * V, Field * ETA, double dt){
	double coriolis = 0;
	double tidalForce = 0;
	double dSurfLat = 0;
	double surfHeight = 0;

	double northEta = 0;
	double southEta = 0;

	double lat, lon;

	for (int i = 0; i < v->fieldLatLen; i++) {
		lat = v->lat[i] * radConv;
		for (int j = 0; j < v->fieldLonLen; j++) {
			lon = v->lon[j] * radConv;

			northEta = ETA->CenterP(i, j);
			southEta = ETA->SouthP(i, j);

			dSurfLat = (northEta - southEta) / (eta->dLat*radConv);

			surfHeight = (consts->g.Value() / consts->radius.Value())*dSurfLat;

			coriolis = 2. * consts->angVel.Value()*sin(lat)*U->SouthWestAvg(i, j);

			tidalForce = consts->loveReduct.Value()*(1. / consts->radius.Value())* dUlat->solution[i][j];

			double val = (-coriolis - surfHeight + tidalForce - V->solution[i][j] * consts->alpha.Value())*dt + VOLD->solution[i][j];

			VNEW->solution[i][j] = (-coriolis - surfHeight + tidalForce - V->solution[i][j] * consts->alpha.Value())*dt + VOLD->solution[i][j];
		}
	}
	
}

void Solver::UpdateSurfaceHeight(Field * ETAOLD, Field * ETANEW, Field * U, Field * V, Field * ETA, double dt){
	double vGrad = 0;
	double uGrad = 0;
	double northv = 0;
	double southv = 0;
	double eastu = 0;
	double westu = 0;

	double lat, lon;


	for (int i = 1; i < eta->fieldLatLen-1; i++) {
		lat = eta->lat[i] * radConv;
		for (int j = 0; j < eta->fieldLonLen; j++) {
			lon = eta->lon[j] * radConv;

			northv = V->CenterP(i - 1, j)*cos(v->lat[i - 1] * radConv);
			southv = V->CenterP(i, j)*cos(v->lat[i] * radConv);

			vGrad = (northv - southv) / (v->dLat*radConv);

			eastu = U->CenterP(i, j);
			westu = U->WestP(i, j);

			uGrad = (eastu - westu) / (u->dLon*radConv);

			ETANEW->solution[i][j] = consts->h.Value() / (consts->radius.Value()*cos(lat))*(-vGrad - uGrad)*dt + ETAOLD->solution[i][j];
		}
	}
	
	for (int j = 0; j < eta->fieldLonLen; j++) {
		ETANEW->solution[0][j] = lagrangeInterp(ETANEW, 0, j);
		ETANEW->solution[eta->fieldLatLen - 1][j] = lagrangeInterp(ETANEW, eta->fieldLatLen - 1, j);
	}

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
	outstring << "Time step: \t" << consts->timeStep.Value() << " seconds\n\n";
	Out->Write(OUT_MESSAGE, &outstring);

	double lat = 0;
	double lon = 0;

	double surfaceheight = 0;
	int output = 0;

	
	int outputTime = 34560;
	int inc = outputTime;
	DumpSolutions(-1);

	//Update mass before calculations begin
	mass->UpdateMass();

	//Update cell energies and globally averaged energy
	energy->UpdateKinE();
	
	while (simulationTime <= consts->endTime.Value() && !energy->converged) {

		UpdatePotential();

		simulationTime = simulationTime + consts->timeStep.Value();
		
		//Solve for v
		UpdateNorthVel(v, vNew, u, v, eta, consts->timeStep.Value());

		//solve for u
		UpdateEastVel(u, uNew, u, v, eta, consts->timeStep.Value());

		//Solve for eta
		UpdateSurfaceHeight(eta, etaNew, uNew, vNew, eta, consts->timeStep.Value());

		//Overwite previous timestep solutions at end of iteration.
		*eta = *etaNew;
		*u = *uNew;
		*v = *vNew;

		if (fmod(iteration, (outputTime / inc)) == 0) {
			//Update mass
			//mass->UpdateMass();

			//Update Kinetic Energy Field
			energy->UpdateKinE();

			//Update Globally Avergaed Kinetic Energy 
			energy->UpdateDtKinEAvg();
		}

		//Check for output
		if (fmod(iteration, outputTime) == 0) {//11520/10/2outputTime/10/2
			energy->UpdateOrbitalKinEAvg(inc);
			//std::cout << energy->isConverged;

			//Check for convergence
			if (orbitNumber>1) energy->IsConverged();

			outstring << std::fixed << std::setprecision(2) << simulationTime / 86400.0 << " days: \t" << 100 * (simulationTime / consts->endTime.Value()) << "%\t" << output;
			Out->Write(OUT_MESSAGE, &outstring);
			
			DumpSolutions(output);
			output++;
		}
		if (fmod(iteration, consts->period.Value()/consts->timeStep.Value()) == 0) orbitNumber++;
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
		CreateDirectory(&(consts->path + "\\EastVelocity\\")[0], NULL);
		CreateDirectory(&(consts->path + "\\NorthVelocity\\")[0], NULL);
		CreateDirectory(&(consts->path + "\\Displacement\\")[0], NULL);
		CreateDirectory(&(consts->path + "\\Grid\\")[0], NULL);
		
		remove(&(consts->path + "\\diss.txt")[0]);

		std::ofstream uLat(consts->path+"\\Grid\\u_lat.txt", std::ofstream::out);
		std::ofstream uLon(consts->path+"\\Grid\\u_lon.txt", std::ofstream::out);
		std::ofstream vLat(consts->path + "\\Grid\\v_lat.txt", std::ofstream::out);
		std::ofstream vLon(consts->path + "\\Grid\\v_lon.txt", std::ofstream::out);

		for (int j = 0; j < u->ReturnFieldLonLen(); j++) uLon << u->lon[j] << '\t';
		for (int i = 0; i < u->ReturnFieldLatLen(); i++) uLat << u->lat[i] << '\t';

		for (int j = 0; j < v->ReturnFieldLonLen(); j++) vLon << v->lon[j] << '\t';
		for (int i = 0; i < v->ReturnFieldLatLen(); i++) vLat << v->lat[i] << '\t';
	}

	else {
		std::string out = std::to_string(out_num);

		std::ofstream uFile(consts->path + "\\EastVelocity\\u_vel_" + out + ".txt", std::ofstream::out);
		std::ofstream vFile(consts->path + "\\NorthVelocity\\v_vel_" + out + ".txt", std::ofstream::out);
		std::ofstream etaFile(consts->path + "\\Displacement\\eta_" + out + ".txt", std::ofstream::out);
		
		//std::ofstream dissFile(consts->path + "diss.txt", std::ofstream::app);

		FILE * dissFile;

		errno_t err = fopen_s(&dissFile, &(consts->path + "\\diss.txt")[0], "a");

		//fprintf(pFile, "Name %d [%-10.10s]\n", n + 1, name);
		fprintf(dissFile, "%.15f \n", energy->orbitDissEAvg[energy->orbitDissEAvg.size() - 1]);

		err = fclose(dissFile);
		//dissFile << energy->orbitDissEAvg[energy->orbitDissEAvg.size() - 1] << std::endl;

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
	}
};

void Solver::ReadInitialConditions(void) {
	std::ifstream eastVel(consts->path + "\\InitialConditions\\u_vel.txt", std::ifstream::in);
	std::ifstream northVel(consts->path + "\\InitialConditions\\v_vel.txt", std::ifstream::in);
	std::ifstream displacement(consts->path + "\\InitialConditions\\eta.txt", std::ifstream::in);

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