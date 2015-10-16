#ifndef SOLVER_H
#define SOLVER_H

#include "field.h"
#include "globals.h"
#include "mesh.h"
#include "energy.h"
#include "outFiles.h"

class Solver {
private:
	int solverType;
	int dumpTime;
	enum Potential {OBLIQ, ECC_RAD, ECC_LIB, ECC, FULL};

	void UpdateEccRadPotential(void);
	void UpdateEccLibPotential(void);
	void UpdateEccPotential(void);
	void UpdateObliqPotential(void);
	void UpdateFullPotential(void);

	void ReadInitialConditions(void);
	void CopyInitialConditions(std::ifstream & file, Field *);

public:
	double simulationTime = 0;
	int orbitNumber = 0;
	int iteration = 0;
	int convergeCount = 0;
	int convergeMax = 5;
	double outputCount = 0.;

	double dt;

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
	Energy * energy;

	Field * etaNew;
	Field * vNew;
	Field * uNew;

	Field * etaNewHalf;
	Field * vNewHalf;
	Field * uNewHalf;

	Field * vDissTerm;
	Field * uDissTerm;

	Field * vNorthEastAvg;
	Field * uSouthWestAvg;

	double * cosMinusB;
	double * cosPlusB;
	double * sinMinusB;
	double * sinPlusB;

	Solver(int type, int dump, Globals *, Mesh *, Field *, Field*, Field *, Field *, Field *, Energy *);

	void Solve();

	int InitialConditions(void);

	void Explicit();
	void Implicit();
	void CrankNicolson();

	void UpdatePotential();
	void UpdateEastVel(Field * U, Field * UNEW, Field * V, Field * ETA);
	void UpdateNorthVel(Field * V, Field * VNEW, Field * U, Field * ETA);
	void UpdateSurfaceHeight(Field * ETA, Field * ETANEW, Field * U, Field * V);

	void InterpPole(Field * Field);

	void DumpSolutions(int output_num);
	void DumpFields(int output_num);
};

#endif
