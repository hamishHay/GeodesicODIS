#ifndef SOLVER_H
#define SOLVER_H

#include "field.h"
#include "globals.h"
#include "mesh.h"
#include "depth.h"
#include "energy.h"
#include "outFiles.h"

#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

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
	Depth * newRadius;
	Energy * energy;

	double ** newRadiusArray;

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

	double ** etaVAvgArray;
	double ** etaUAvgArray;
	double ** etaInterpArray;

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

	Field * etaVAvg;
	Field * etaUAvg;
	Field * etaInterp;

	double * cosMinusB;
	double * cosPlusB;
	double * sinMinusB;
	double * sinPlusB;

	double ** SH_cos_coeff;
	double ** SH_sin_coeff;

	int l_max;

  //------------------------Objects for HDF5 Storage----------------------------

  float * eta_1D;
  float * u_1D;
  float * v_1D;

  hsize_t eta_rows;
  hsize_t eta_cols;
  hsize_t u_rows;
  hsize_t u_cols;
  hsize_t v_rows;
  hsize_t v_cols;

  static hsize_t rank_field = 3;

  const char dataFilePath[] = "data.h5";

  hsize_t start[3];
  hsize_t count[3];

  hid_t file;

  hsize_t dims[rank_field];
  hsize_t max_dims[rank_field];

  hid_t data_space;
  hid_t mem_space;
  hid_t data_set;

  //
  // const H5std_string FILE_NAME("data.h5"); // File name for all solvable variables
  //
  // const
  //
  // const H5std_string DATASET_NAME("displacement");
  // const H5std_string DATASET_NAME("north velocity");
  // const H5std_string DATASET_NAME("east velocity");
  //
  // float * etaOutput;
  // float * vOutput;
  // float * uOutput;
  //
  // H5File ODIS_Data( FILE_NAME, H5F_ACC_TRUNC); // Master data file
  //
  // hsize_t dimsf_eta[2] = {etaLatLen, etaLonLen};
  // DataSpace dataspace_eta(2, dimsf_eta);

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
	inline void InterpSurfaceHeight() __attribute__((always_inline));

	inline void ExtractSHCoeff() __attribute__((always_inline));

	void InterpPole(Field * Field);

	void DumpSolutions(int output_num, double time);
	void DumpFields(int output_num);
};

#endif
