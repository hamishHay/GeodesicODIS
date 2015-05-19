#include "solver.h"
#include "globals.h"
#include "field.h"
#include "mesh.h"
#include "mathRoutines.h"
#include "mass.h"
#include "energy.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdio>

Solver::Solver(int type, int dump, Globals * Consts, Mesh * Grid, Field * UGradLon, Field * UGradLat, Field * VelU, Field * VelV, Field * Eta, Mass * MassField, Energy * EnergyField) {
	solverType = type;
	dumpTime = dump;

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

};

void Solver::InitialConditions(int action) {
	if (action) {
		UpdateEccPotential();
	}
	else {
		//Read file
	}

}

void Solver::Solve() {
	std::cout << "\nEntering solver: ";
	switch (solverType) {
	case 0:
		std::cout << "Explicit\n";
		Explicit();
		break;
	default:
		std::cout << "No solver type selected.\n Terminating execution.";
	}
};

void Solver::UpdateEccPotential() {
	double lat = 0;
	double lon = 0;
	double B = consts->angVel.Value()*simulationTime;

	//CHECK FOR ECCENTRICITY OR OBLIQUITY!
	double potential = 0.25*pow(consts->angVel.Value(), 2)*pow(consts->radius.Value(), 2)*consts->theta.Value();

	for (int i = 0; i < dUlat->fieldLatLen; i++) {
		lat = dUlat->lat[i] * radConv;
		for (int j = 0; j < dUlat->fieldLonLen; j++) {
			lon = dUlat->lon[j] * radConv;
			//dUlat->solution[i][j] = potential*(-9.*sin(2*lat)*cos(B)); // P20
			//dUlat->solution[i][j] = potential*(- 1.5*sin(2*lat) * (7 * cos(2*lon - B) - cos(2*lon + B))); //P22 Negative sign necessary because of diff between finite diff grad sign and analytic sign
			dUlat->solution[i][j] = -6 * potential*cos(2 * lat)*(cos(lon - B) + cos(lon + B)); //P21 cos(2*lat)*
			//dUlat->solution[i][j] = 0.25*pow(consts->angVel, 2)*pow(consts->radius, 2)*consts->theta*-6 *cos(2 * lat)*(cos(lon - B) + cos(lon + B));
			//dUlat->solution[i][j] += 0.25*pow(consts->angVel, 2)*pow(consts->radius, 2)*consts->e*(-9.*sin(2 * lat)*cos(B));
			//dUlat->solution[i][j] += 0.25*pow(consts->angVel, 2)*pow(consts->radius, 2)*consts->e*(-1.5*sin(2 * lat) * (7 * cos(2 * lon - B) - cos(2 * lon + B)));
		}
	}
		// Skip for P20 only as dUlon = 0;
		for (int i = 0; i < dUlon->fieldLatLen; i++) {
			lat = dUlon->lat[i] * radConv;
			for (int j = 0; j < dUlon->fieldLonLen; j++) {
				lon = dUlon->lon[j] * radConv;
				//dUlon->solution[i][j] = potential *3* pow(cos(lat), 2)*(-7 * sin(2*lon - B) + sin(2*lon + B)); //P22
				dUlon->solution[i][j] = 3 * potential * sin(2*lat)*(sin(lon - B) + sin(lon + B)); //P21
				//dUlon->solution[i][j] = 0.25*pow(consts->angVel, 2)*pow(consts->radius, 2)*consts->theta * 3 * sin(2 * lat)*(sin(lon - B) + sin(lon + B));
				//dUlon->solution[i][j] += 0.25*pow(consts->angVel, 2)*pow(consts->radius, 2)*consts->e * 3 * pow(cos(lat), 2)*(-7 * sin(2 * lon - B) + sin(2 * lon + B));
			}
		}
}

void Solver::UpdateEastVel(int i, int j, double lat, double lon){
	double coriolis = 0;
	double tidalForce = 0;

	double dSurfLon = 0;
	double surfHeight = 0;

	double eastEta = 0;
	double westEta = 0;

	//eastEta = (eta->EastP(i, j) + eta->EastP(i, j + 1))*0.5;
	//westEta = (eta->CenterP(i, j) + eta->WestP(i, j))*0.5;
	eastEta = eta->EastP(i, j);
	westEta = eta->CenterP(i, j);


	dSurfLon = (eastEta - westEta) / (eta->dLon*radConv);

	surfHeight = (consts->g.Value() / (consts->radius.Value()*cos(lat)))*dSurfLon;

	coriolis = 2. * consts->angVel.Value()*sin(lat)*v->NorthEastAvg(i, j);
	tidalForce = consts->loveReduct.Value()*(1. / (consts->radius.Value()*cos(lat)))* dUlon->solution[i][j];
	uNew->solution[i][j] = (coriolis - surfHeight + tidalForce - u->solution[i][j] * consts->alpha.Value())*consts->timeStep.Value() + u->solution[i][j];
}

void Solver::UpdateNorthVel(int i, int j, double lat, double lon){
	double coriolis = 0;
	double tidalForce = 0;
	double dSurfLat = 0;
	double surfHeight = 0;

	double northEta = 0;
	double southEta = 0;
	
	/*
	if (i == 0){
		northEta = (eta->CenterP(i, j) + eta->CenterP(i, eta->opp[j]))*0.5;
		//northEta = eta->CenterP(i, j);
		//southEta = eta->SouthP(i, j);
		//dSurfLat = (northEta - southEta) / (eta->dLat*radConv);

	}
	else {
		northEta = (eta->CenterP(i, j) + eta->NorthP(i, j))*0.5;
	}
	if (i == v->fieldLatLen - 1){
		southEta = (eta->CenterP(i + 1, j) + eta->SouthP(i+1, j))*0.5;
		//northEta = eta->CenterP(i,j)
	}
	else southEta = (eta->SouthP(i, j) + eta->SouthP(i + 1, j))*0.5;*/

	northEta = eta->CenterP(i, j);
	southEta = eta->SouthP(i, j);

	dSurfLat = (northEta - southEta) / (eta->dLat*radConv);
	
	surfHeight = (consts->g.Value() / consts->radius.Value())*dSurfLat;
	if (i == v->fieldLatLen - 1) {
		coriolis = 2. * consts->angVel.Value()*sin(lat)*(u->CenterP(i, j) + u->WestP(i, j))*0.5;
	}
	else coriolis = 2. * consts->angVel.Value()*sin(lat)*u->SouthWestAvg(i, j);
	
	tidalForce = consts->loveReduct.Value()*(1. / consts->radius.Value())* dUlat->solution[i][j];
	vNew->solution[i][j] = (-coriolis - surfHeight + tidalForce - v->solution[i][j] * consts->alpha.Value())*consts->timeStep.Value() + v->solution[i][j];
}

void Solver::UpdateSurfaceHeight(int i, int j, double lat, double lon){
	double vGrad = 0;
	double uGrad = 0;
	double northv = 0;
	double southv = 0;
	double eastu = 0;
	double westu = 0;

	if (i == 0) {
		northv = vNew->NorthP(i, j)*cos(v->lat[i] * radConv);
		//northv = (vNew->NorthP(i, j)*cos(v->lat[i] * radConv) + vNew->NorthP(i - 1, j)*cos(v->lat[i+1] * radConv))*0.5;
	}
	//else if (i == 1) {
	//	northv = (vNew->CenterP(i - 1, j)*cos(v->lat[i - 1] * radConv) + vNew->NorthP(i - 1, j)*cos(v->lat[i - 1] * radConv))*0.5;
	//}
	else {
		northv = vNew->CenterP(i-1, j)*cos(v->lat[i - 1] * radConv);
		//northv = (vNew->CenterP(i - 1, j)*cos(v->lat[i - 1] * radConv) + vNew->NorthP(i - 1, j)*cos(v->lat[i - 2] * radConv))*0.5;
	}

	if (i == v->fieldLatLen - 1) { 
		southv = vNew->CenterP(i, j)*cos(v->lat[i] * radConv);
		//southv = (vNew->CenterP(i, j)*cos(v->lat[i] * radConv) + vNew->SouthP(i, j)*cos(v->lat[i] * radConv))*0.5;
	}
	else if (i == v->fieldLatLen) {
		southv = vNew->SouthP(i -1 , j)*cos(v->lat[i - 1] * radConv);
		//southv = (vNew->SouthP(i-1, j)*cos(v->lat[i - 1] * radConv) + vNew->SouthP(i-2, j)*cos(v->lat[i - 2] * radConv))*0.5;
	}
	else {
		southv = vNew->CenterP(i, j)*cos(v->lat[i] * radConv);
		//southv = (vNew->CenterP(i, j)*cos(v->lat[i] * radConv)+vNew->SouthP(i,j)*cos(v->lat[i+1]*radConv))*0.5;
	}

	vGrad = (northv - southv)/(v->dLat*radConv);

	//eastu = (uNew->CenterP(i, j) + uNew->EastP(i, j))*0.5;
	//westu = (uNew->WestP(i, j) + uNew->WestP(i, j-1))*0.5;

	eastu = uNew->CenterP(i, j);
	westu = uNew->WestP(i, j);

	uGrad = (eastu - westu) / (u->dLon*radConv);

	etaNew->solution[i][j] = consts->h.Value() / (consts->radius.Value()*cos(lat))*(-vGrad - uGrad)*consts->timeStep.Value() + eta->solution[i][j];

}


void Solver::Explicit() {
	//Check for stability
	std::cout << "Entering time loop:\n\n";
	std::cout << "End time: \t" << consts->endTime.Value() / 86400.0 << " days\n";
	std::cout << "Time step: \t" << consts->timeStep.Value() << " seconds\n\n";

	double lat = 0;
	double lon = 0;

	double surfaceheight = 0;
	int output = 0;

	
	int outputTime = 46080;
	int inc = outputTime;
	DumpSolutions(-1);

	//Update mass before calculations begin
	mass->UpdateMass();

	//Update cell energies and globally averaged energy
	energy->UpdateKinE();
	
	while (simulationTime <= consts->endTime.Value()) {

		UpdateEccPotential();

		simulationTime = simulationTime + consts->timeStep.Value();
		
		//Solve for v
		for (int i = 1; i < v->fieldLatLen-1; i++) {
			lat = v->lat[i]* radConv;
			for (int j = 0; j < v->fieldLonLen; j++) {
				lon = v->lon[j] * radConv;
				UpdateNorthVel(i, j, lat, lon);
			}
		}
		
		for (int j = 0; j < v->fieldLonLen; j++) {
			vNew->solution[0][j] = lagrangeInterp(vNew, 0, j);
			vNew->solution[v->fieldLatLen - 1][j] = lagrangeInterp(vNew, v->fieldLatLen - 1, j);
			//vNew->solution[0][j] = linearInterp(vNew, 0, j);
			//vNew->solution[v->fieldLatLen - 1][j] = linearInterp(vNew, v->fieldLatLen - 1, j);
		}
		
		//solve for u
		for (int i = 1; i < u->fieldLatLen-1; i++) {
			lat = u->lat[i]* radConv;
			for (int j = 0; j < u->fieldLonLen; j++) {
				lon = u->lon[j] * radConv;
				UpdateEastVel(i, j, lat, lon);
			}
		}
		/*
		for (int j = 0; j < eta->fieldLonLen; j++) {
			//uNew->solution[1][j] = lagrangeInterp(uNew, 1, j);
			//uNew->solution[u->fieldLatLen - 2][j] = lagrangeInterp(uNew, u->fieldLatLen - 2, j);
			//vNew->solution[0][j] = linearInterp(vNew, 0, j);
			//vNew->solution[v->fieldLatLen - 1][j] = linearInterp(vNew, v->fieldLatLen - 1, j);
		}*/

		
		//Solve for eta 
		for (int i = 1; i < eta->fieldLatLen-1; i++) {
			lat = eta->lat[i] * radConv;
			for (int j = 0; j < eta->fieldLonLen; j++) {
				lon = eta->lon[j] * radConv;
				UpdateSurfaceHeight(i, j, lat, lon);
			}
		}
		
		for (int j = 0; j < eta->fieldLonLen; j++) {
			etaNew->solution[0][j] = lagrangeInterp(etaNew, 0, j);
			etaNew->solution[eta->fieldLatLen - 1][j] = lagrangeInterp(etaNew, eta->fieldLatLen - 1, j);
			//etaNew->solution[0][j] = linearInterp(etaNew, 0, j);
			//etaNew->solution[eta->fieldLatLen - 1][j] = linearInterp(etaNew, eta->fieldLatLen - 1, j);
		}
		
		//Average eta out at poles
		
		double npoleSum = 0;
		double spoleSum = 0;
		for (int j = 0; j < eta->fieldLonLen; j++) {
			npoleSum += etaNew->solution[0][j];
			spoleSum += etaNew->solution[etaNew->fieldLatLen - 1][j];
		}
		npoleSum = npoleSum / etaNew->fieldLonLen;
		spoleSum = spoleSum / etaNew->fieldLonLen;

		for (int j = 0; j < eta->fieldLonLen; j++) {
			etaNew->solution[0][j] = npoleSum;
			etaNew->solution[etaNew->fieldLatLen - 1][j] = spoleSum;
		}

		//Overwite previous timestep solutions at end of iteration.
		*eta = *etaNew;
		*u = *uNew;
		*v = *vNew;

		

		if (fmod(iteration, (outputTime / inc)) == 0) {
			//Update mass
			mass->UpdateMass();

			//Update Kinetic Energy Field
			energy->UpdateKinE();

			//Update Globally Avergaed Kinetic Energy 
			energy->UpdateDtKinEAvg();
		}

		//Check for output
		if (fmod(iteration, outputTime) == 0) {//11520/10/2
			energy->UpdateOrbitalKinEAvg(inc);

			std::cout << std::fixed << std::setprecision(2) << simulationTime / 86400.0 << " days: \t" << 100 * (simulationTime / consts->endTime.Value()) << "%\t" << output << std::endl;
			DumpSolutions(output);
			output++;
		}
		if (simulationTime >= iteration*consts->period.Value()) orbitNumber++;
		iteration++;
	}

	

	//deallocate memory assigned to temporary fields
	delete etaNew;
	delete uNew;
	delete vNew;

	std::cout << "\nSimulation complete.\n\nEnd time: \t\t" << simulationTime / 86400.0 << "\n";
	std::cout << "Total iterations: \t" << iteration;
};

void Solver::DumpSolutions(int out_num) {

	std::string path = "C:\\Users\\Hamish\\Google Drive\\LPL\\Icy Satellites\\Results\\Working\\";
	if (out_num == -1) {

		remove("C:\\Users\\Hamish\\Google Drive\\LPL\\Icy Satellites\\Results\\Working\\diss.txt");

		std::ofstream uLat(path+"u_lat.txt", std::ofstream::out);
		std::ofstream uLon(path+"u_lon.txt", std::ofstream::out);
		std::ofstream vLat(path + "v_lat.txt", std::ofstream::out);
		std::ofstream vLon(path + "v_lon.txt", std::ofstream::out);

		for (int j = 0; j < u->ReturnFieldLonLen(); j++) uLon << u->lon[j] << '\t';
		for (int i = 0; i < u->ReturnFieldLatLen(); i++) uLat << u->lat[i] << '\t';

		for (int j = 0; j < v->ReturnFieldLonLen(); j++) vLon << v->lon[j] << '\t';
		for (int i = 0; i < v->ReturnFieldLatLen(); i++) vLat << v->lat[i] << '\t';
	}

	else {
		std::string out = std::to_string(out_num);

		//std::string uVelString = "u_vel_" + out + ".txt";

		std::ofstream uFile(path + "u_vel_" + out + ".txt", std::ofstream::out);
		std::ofstream vFile(path + "v_vel_" + out + ".txt", std::ofstream::out);
		//std::ofstream uLat("u_lat.txt", std::ofstream::out);
		//std::ofstream uLon("u_lon.txt", std::ofstream::out);
		std::ofstream etaFile(path + "eta_" + out + ".txt", std::ofstream::out);
		//std::ofstream dULonFile("dU_lon_" + out + ".txt", std::ofstream::out);

		
		std::ofstream dissFile(path + "diss.txt", std::ofstream::app);

		dissFile << energy->orbitDissEAvg[energy->orbitDissEAvg.size() - 1] << std::endl;

		//for (int j = 0; j < u->ReturnFieldLonLen(); j++) {
		//	uLon << u->grid->lon[j * 2] << '\t';
		//}

		for (int i = 0; i < u->ReturnFieldLatLen(); i++) {
			//uLat << u->grid->lat[i * 2] << '\t';
			for (int j = 0; j < u->ReturnFieldLonLen(); j++) {
				uFile << uNew->solution[i][j] << '\t';
				etaFile << etaNew->solution[i][j] << '\t';
			}
			uFile << std::endl;
			etaFile << std::endl;
		}

		for (int i = 0; i < v->fieldLatLen; i++) {
			//uLat << u->grid->lat[i * 2] << '\t';
			for (int j = 0; j < v->fieldLonLen; j++) {
				vFile << vNew->solution[i][j] << '\t';
			}
			vFile << std::endl;
		}

		
	}
}

