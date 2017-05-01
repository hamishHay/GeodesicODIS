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

/* Solver class
 *
 * This class is designed to contain the methods by which the Laplace Tidal
 * Equations can be solved. Currently the only solver implemented is a simple
 * first-order explicit Euler time-scheme.
 *
 * The solver contains pointers to all arrays relevant to the solution of the
 * Laplace Tidal Equations, including velocity components, tidal displacement,
 * and forcings. It also contains the relevant functions for saving data into
 * compressed hdf5 files.
 *
 * For further detailed documention of the class methods see "solver.cpp".
*/

class Solver {
private:
  int solverType;                // Type of solver selected in main
  int dumpTime;                  // TODO - REMOVE UNUSED PARAM

  // enum for different types of tidal forcing
  enum Potential {OBLIQ,         // Obliquity tide (Tyler 2011)
                  ECC_RAD,       // Eccentricity radial tide (Tyler 2011)
                  ECC_LIB,       // Eccentricity libration tide (Tyler 2011)
                  ECC,           // ECC_RAD + ECC_LIB
                  FULL,          // ECC_RAD + ECC_LIB + OBLIQ
                  TOTAL,         // Entire ecc and obliq potential to second order in Eccentricity and Obliquity
                  ECC_W3,        // Time-dependent degree-3 eccentricity tide // TODO - Add expressions to documentation
                  OBLIQ_W3};     // Time-dependent degree-3 obliquity tide

  void ReadInitialConditions(bool); // Reads hdf5 init condition file for u,v and eta.

public:
  double simulationTime = 0;              // Current time in simulation
  int orbitNumber = 0;                    // Current orbit in simulation
  int iteration = 0;
  int convergeCount = 0;                  // Number of full orbits since convergence criteria met
  int convergeMax = 5;                    // Number of orbits to continue simulation since convergence
  int output = 0;                         // Number of current data output
  double outputCount = 0.;                // Simulation time between data output events
  bool loading;                           // Switch for using ocean loading and self-attraction (Matsuyama 2014)
  bool bc;

  double dt;                              // Time-step

  std::ostringstream outstring;           // String for writing output messages

  OutFiles * Out;                         // Pointer to all verbose output and simulation path

  Potential tide;                         // Enum for tidal forcing, see above.

  Globals * consts;                       // Pointer to general simulation constants
  Mesh * grid;                            // Pointer to full grid information
  Field * dUlon;                          // Pointer to field potential gradient in longitude field
  Field * dUlat;                          // Pointer to field potential gradient in latitude field
  Field * eta;                            // Pointer to displacement solution
  Field * v;                              // Pointer to north velocity solution
  Field * u;                              // Pointer to east velocity solution
  Depth * depth;                          // Pointer to ocean depth constant
  Energy * energy;                        // Pointer to kinetic energy solution

  double ** cellArea;

  /* Each solvable field requires information to be stored at the new and
   * current timesteps. These are denotated as Old (current dt) and New (next dt).
   * "Array" is appended to each of the following double pointers to distinguish
   * them from the equivalent Field instances. Solver accesses the 2D arrays
   * directly, rather than through their classes, to avoid overhead with
   * accessing class member variables and member functions.
  */

  double ** etaOldArray;
  double ** etaNewArray;

  double ** vOldArray;
  double ** vNewArray;

  double ** uOldArray;
  double ** uNewArray;

  double ** dUlonArray;
  double ** dUlatArray;

  double ** depthArray;                   // Array of ocean thickness

  double ** oceanLoadingArrayU;           // Derivative of the ocean loading term at u nodes
  double ** oceanLoadingArrayV;           // Derivative of the ocean loading term at v nodes
  double ** oceanLoadingArrayEta;

  double ** vDissArray;                   // Array to contain v dissipation term across the gird
  double ** uDissArray;                   // Array to contain u dissipation term across the gird


  // The numerical grid in ODIS is staggerd in space. Velocity solutions thus
  // require averages of each other to determine the solutions at each timestep.
  // The averages are used in both the coriolis terms, and the bottom drag terms
  // if using. (see Hay and Matsuyama, 2017).
  double ** vNEAvgArray;                  // North East average of v (req'd to solve for u)
  double ** uSWAvgArray;                  // South West average of u (req'd to solver for v)

  double ** etaUAvgArray;
  double ** etaVAvgArray;



  /* To include the effects of ocean loading and self-attraction (Matsuyama 2014),
   * it is necessary to compute the derivative of the displacement field when it
   * is expanded into spherical harmonic components. These derivatives are used
   * to update the velocity fields. For the north velocity component, derivatives
   * of the associated legendre polynomials are required. These calculated and
   * stored for every latitude, degree l and order m in the 3D pointer
   * vdLegendreArray. Derivatives in longitude require only the associated
   * Legendre function, which is stored in uLegendreArray at every latitude for
   * degree l and order m.
  */
  double *** vdLegendreArray;             // Derivative of associated legendre polynomials at v nodes
  double *** uLegendreArray;              // Associated legendre polynomials at u nodes
  // double *** etaLegendreArray;
  double * etaLegendreArray;


  double ** vCosMLon;                     // cos(m*lon) at v nodes ([m][lon])
  double ** vSinMLon;                     // sin(m*lon) at v nodes ([m][lon])
  double ** uCosMLon;                     // cos(m*lon) at u nodes ([m][lon])
  double ** uSinMLon;                     // sin(m*lon) at u nodes ([m][lon])
  double * etaCosMLon;                   // cos(m*lon) at eta nodes ([m][lon])
  double * etaSinMLon;                   // sin(m*lon) at eta nodes ([m][lon])

  double ** boundaryArray;


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

  double radius;
  double angVel;
  double theta;
  double ecc;
  double smAxis;

  // THESE ARRAYS ARE INEFFICIENT AND SO SHOULD BE REMOVED TODO: DELETE
  double * cosMinusB;
  double * cosPlusB;
  double * sinMinusB;
  double * sinPlusB;

  double * SH_cos_coeff;
  double * SH_sin_coeff;

  double * loadK;
  double * loadH;
  double * gammaFactor;

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
  int UpdateLoading();
  int ApplyBoundaries();
  inline void InterpSurfaceHeight() __attribute__((always_inline));

  int ExtractSHCoeff();
  int LegendreDeriv();
  int Legendre();
  int InterpPoles();
  int FindAverages();

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
