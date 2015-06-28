#ifndef SOLVER_H
#define SOLVER_H

#include "field.h"
#include "globals.h"
#include "mesh.h"
#include "mass.h"
#include "energy.h"
#include "outFiles.h"

class Solver {
private:
	int solverType;
	int dumpTime;
	enum Potential {OBLIQ, ECC_RAD, ECC_LIB, FULL};

	void UpdateEccRadPotential(void);
	void UpdateEccLibPotential(void);
	void UpdateObliqPotential(void);
	void UpdateFullPotential(void);

	void ReadInitialConditions(void);

public:
	double simulationTime = 0;
	int orbitNumber = 1;
	int iteration = 0;

	std::ostringstream outstring;

	OutFiles * Out;

	Potential tide;

	Globals * consts;
	Mesh * grid;
	Field * dUlon;
	Field * dUlat;
	Field * eta;
	Field * v;
	Field * u;
	Mass * mass;
	Energy * energy;

	Field * etaNew;
	Field * vNew;
	Field * uNew;

	Solver(int type, int dump, Globals *, Mesh *, Field *, Field*, Field *, Field *, Field *, Mass *, Energy *);

	void Solve();
	
	int InitialConditions(void);
	
	void Explicit();
	void Implicit();
	void CrankNicolson();

	void UpdatePotential();
	void UpdateEastVel(int i, int j, double lat, double lon);
	void UpdateNorthVel(int i, int j, double lat, double lon);
	void UpdateSurfaceHeight(int i, int j, double lat, double lon);

	void DumpSolutions(int output_num);
};

#endif