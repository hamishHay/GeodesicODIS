#ifndef SOLVER_H
#define SOLVER_H

#include "field.h"
#include "globals.h"
#include "mesh.h"
#include "depth.h"
#include "energy.h"
#include "outFiles.h"

#include "H5Cpp.h"

#include <signal.h>

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

class Solver {
private:
	int solverType;
	int dumpTime;
	enum Potential {OBLIQ, ECC_RAD, ECC_LIB, ECC, FULL, TOTAL, ECC_W3};

	void UpdateEccRadPotential(void);
	void UpdateEccLibPotential(void);
	void UpdateEccPotential(void);
	inline void UpdateObliqPotential(void) __attribute__((always_inline));
	void UpdateFullPotential(void);
	void UpdateTotalPotential(void);
  void UpdateEccDeg3Potential(void);

	void ReadInitialConditions(void);

public:
	double simulationTime = 0;
	int orbitNumber = 0;
	int iteration = 0;
	int convergeCount = 0;
	int convergeMax = 5;
	int output = 0;
	double outputCount = 0.;
  bool loading;

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
  double *** etaLegendreArray;
  double ** etaCosMLon;
  double ** etaSinMLon;

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

  double * loadK;
  double * loadH;
  double * gammaFactor;
  double ** shPower;
  int ** lm_solve;
  int l_solve_len;
  int * m_solve_len;

	int l_max;

	Solver(int type, int dump, Globals *, Mesh *, Field *, Field*, Field *, Field *, Field *, Energy *, Depth *);

	void Solve();

	int InitialConditions(void);

	int Explicit();
	void Implicit();
	void CrankNicolson();

	inline void UpdatePotential() __attribute__((always_inline));
	void UpdateEastVel();
	int UpdateNorthVel();
	void UpdateSurfaceHeight();
	inline void InterpSurfaceHeight() __attribute__((always_inline));

	int ExtractSHCoeff();

	void InterpPole(Field * Field);

	void DumpSolutions(int output_num, double time);
	void DumpFields(int output_num);

  void CreateHDF5FrameWork(void);

  // void CatchExit(int sig, int output);

  //------------------------Objects for HDF5 Storage----------------------------

  float * eta_1D;
  float * u_1D;
  float * v_1D;
  float * diss_1D;
  float * diss_avg_1D;
  float * harm_coeff_1D;

  hsize_t eta_rows;
  hsize_t eta_cols;
  hsize_t u_rows;
  hsize_t u_cols;
  hsize_t v_rows;
  hsize_t v_cols;
  hsize_t harm_rows;
  hsize_t harm_cols;

  hsize_t rank_field;
  hsize_t rank_harm;
  hsize_t rank_1D;

  hid_t file;

  hid_t data_space_eta;
  hid_t data_space_u;
  hid_t data_space_v;
  hid_t data_space_diss;
  hid_t data_space_1D_avg;
  hid_t data_space_harm_coeff;

  hid_t mem_space_eta;
  hid_t mem_space_u;
  hid_t mem_space_v;
  hid_t mem_space_diss;
  hid_t mem_space_1D_avg;
  hid_t mem_space_harm_coeff;

  hid_t data_set_eta;
  hid_t data_set_u;
  hid_t data_set_v;
  hid_t data_set_diss;
  hid_t data_set_1D_avg;
  hid_t data_set_harm_coeff;

  char * dataFilePath;

  hsize_t * dims_eta;
  hsize_t * dims_u;
  hsize_t * dims_v;
  hsize_t * dims_1D_avg;
  hsize_t * dims_harm_coeff;

  hsize_t * max_dims_eta;
  hsize_t * max_dims_u;
  hsize_t * max_dims_v;
  hsize_t * max_dims_1D_avg;
  hsize_t * max_dims_harm_coeff;

  hsize_t * start;
  hsize_t * count;

  hsize_t * start_1D;
  hsize_t * count_1D;

  hsize_t * start_harm;
  hsize_t * count_harm;

  hsize_t time_slices;
};

#endif
