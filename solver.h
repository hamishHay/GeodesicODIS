#ifndef SOLVER_H
#define SOLVER_H

#include "field.h"
#include "globals.h"
#include "mesh.h"
#include "depth.h"
#include "energy.h"
#include "outFiles.h"

class Solver {
private:
	int solverType;
	int dumpTime;
	enum Potential {OBLIQ, ECC_RAD, ECC_LIB, ECC, FULL, TOTAL};

	void UpdateEccRadPotential(void);
	void UpdateEccLibPotential(void);
	void UpdateEccPotential(void);
	inline void UpdateObliqPotential(void) __attribute__((always_inline));
	void UpdateFullPotential(void);
	void UpdateTotalPotential(void);

	void ReadInitialConditions(void);
	void CopyInitialConditions(std::ifstream & file, Field *);

public:
	double simulationTime = 0;
	int orbitNumber = 0;
	int iteration = 0;
	int convergeCount = 0;
	int convergeMax = 5;
	int output = 0;
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
	Depth * depth;
	Energy * energy;

	double ** etaOldArray;
	double ** etaNewArray;

	double ** vOldArray;
	double ** vNewArray;

	double ** uOldArray;
	double ** uNewArray;

	double ** dUlonArray;
	double ** dUlatArray;

	double ** vDissArray;
	double ** uDissArray;

	double ** vNEAvgArray;
	double ** uSWAvgArray;

	double ** depthArray;

	int uLatLen;
	int uLonLen;
	double udLon;
	double udLat;

	int vLatLen;
	int vLonLen;
	double vdLon;
	double vdLat;

	int etaLatLen;
	int etaLonLen;
	double etadLon;
	double etadLat;

	Field * etaOld;
	Field * etaNew;
	Field * vOld;
	Field * vNew;
	Field * uOld;
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

	Solver(int type, int dump, Globals *, Mesh *, Field *, Field*, Field *, Field *, Field *, Energy *, Depth *);

	void Solve();

	int InitialConditions(void);

	void Explicit();
	void Implicit();
	void CrankNicolson();

	inline void UpdatePotential() __attribute__((always_inline));
	inline void UpdateEastVel() __attribute__((always_inline));
	inline void UpdateNorthVel() __attribute__((always_inline));
	inline void UpdateSurfaceHeight() __attribute__((always_inline));

	void InterpPole(Field * Field);

	void DumpSolutions(int output_num);
	void DumpFields(int output_num);
};

#endif
