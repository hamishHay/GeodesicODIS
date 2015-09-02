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
	Energy * energy;

	Field * etaNew;
	Field * vNew;
	Field * uNew;

	Field * etaNewHalf;
	Field * vNewHalf;
	Field * uNewHalf;

	Field * vDissTerm;
	Field * uDissTerm;

	std::vector <double> cosMinusB;
	std::vector <double> cosPlusB;
	std::vector <double> sinMinusB;
	std::vector <double> sinPlusB;

	Solver(int type, int dump, Globals *, Mesh *, Field *, Field*, Field *, Field *, Field *, Energy *);

	void Solve();

	int InitialConditions(void);

	void Explicit();
	void Implicit();
	void CrankNicolson();

	void UpdatePotential();
	void UpdateEastVel(Field * UOLD, Field * UNEW, Field * U, Field * V, Field * ETA, double dt);
	void UpdateNorthVel(Field * VOLD, Field * VNEW, Field * U, Field * V, Field * ETA, double dt);
	void UpdateSurfaceHeight(Field * ETAOLD, Field * ETANEW, Field * U, Field * V, Field * ETA, double dt);

	void InterpPole(Field * Field);

	void DumpSolutions(int output_num);
	void DumpFields(int output_num);
};

#endif
