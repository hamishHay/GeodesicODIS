// FILE: globals.h
// DESCRIPTION: Header file for the globals class. Globals is an object designed
//              to store useful model constants that, if necessary, are accessible
//              anywhere in the code. These are listed and defined below in the
//              public section of the class below, but include useful quantities
//              like planetary radius, rotation rate, drag coefficient, etc. Each
//              variable is of type GlobalVar (defined in globalVar.h), with the
//              corresponding machine type held therein.
//
//   Version              Date                   Programmer
//     0.1      |        13/09/2016        |        H. Hay        |
// ---- Initial version of ODIS, described in Hay and Matsuyama (2016)


#ifndef GLOBALS_H
#define GLOBALS_H

#ifdef _WIN32
#define SEP "\\"

#elif _WIN64
#define SEP "\\"

#elif __linux__
#define SEP "/"

#else
#error "OS not supported!"
#endif

#include <string>
#include <array>
#include <vector>
#include "globalVar.h"
#include "outFiles.h"

class OutFiles;           // Forward declare OutFiles class

// define pi here - c++ has no built in value of pi, so it is explicitly defined
// here
const double pi = 3.1415926535897932384626433832795028841971693993751058;
const double radConv = pi / 180.0; // global constant for converting deg --> rad
const int PATH = 1028; // global int for number of characters in the system path

enum Friction { LINEAR,
                QUADRATIC }; // enumerate for drag type, selected by user

enum Surface { FREE,        // Free surface --> no ice lid
               LID,         // Lid surface --> moveable lid
               INF_LID};    // Lid surface --> imoveable (infite rigidity) lid

enum Solver { EULER,        // Explicit Euler time integration
              AB3,          // Adams-Bashforth 3rd order explicit time integration
              RK4 };        // Runge-Kutta 4th order time integration

enum Potential {OBLIQ,      // Obliquity tide (Tyler 2011)
                OBLIQ_WEST, // Westward component of obliquity tide (Tyler 2011)
                OBLIQ_EAST, // Eastward component of obliquity tide (Tyler 2011)
                ECC_RAD,    // Eccentricity radial tide (Tyler 2011)
                ECC_LIB,    // Eccentricity libration tide (Tyler 2011)
                ECC,        // ECC_RAD + ECC_LIB
                ECC_WEST,
                ECC_EAST,
                FULL,       // ECC_RAD + ECC_LIB + OBLIQ
                TOTAL,      // Entire ecc and obliq potential to second order in Eccentricity and Obliquity
                ECC_W3,     // Time-dependent degree-3 eccentricity tide // TODO - Add expressions to documentation
                OBLIQ_W3};  // Time-dependent degree-3 obliquity tide


//Class stores all global values
class Globals {
private:
  // Two member functions for setting all constants (default is for Titan)
  void SetDefault(void);
  void OutputConsts(void); // Print all constants to output file.

public:
  // Class containing functions to output errors and messages, or terminate ODIS.

  OutFiles * Output;

  // Vector to store a 'list' of all global constants.
  // Note that it stores the template IGlobalVar pointer rather than GlobalVar.
  // This is allows multiple types to be stored in one vector instance.
  std::vector<IGlobalVar *> allGlobals;

  Friction fric_type;
  Solver solver_type;
  Potential tide_type;
  Surface surface_type;

  // Stringstream for composing outgoing messages, which is then passed to Output.
  std::ostringstream outstring;

  // string to store path of current simulation.
  std::string path;

  // string for each output tag
  std::vector<std::string> out_tags;

  // character array for identical string from 'path'
  char cpath[PATH];

  int node_num;

  GlobalVar<double> angVel; // Angular velocity
  GlobalVar<double> radius; // Body radius
  GlobalVar<double> loveK2; // k_2 Love Number
  GlobalVar<double> loveH2; // h_2 Love Number
  GlobalVar<double> loveReduct; //Love's reduction factor
  GlobalVar<double> h; //Ocean thickness
  GlobalVar<double> g; //Surface Gravity
  GlobalVar<double> a; //SemiMajor Axis
  GlobalVar<double> e; //Eccentricity
  GlobalVar<double> theta; // obliquity
  GlobalVar<double> timeStep; // Simulation timestep
  GlobalVar<double> alpha; // drag coefficient
  GlobalVar<int> l_max; //Maximum spherical harmonic degree for expansions
  GlobalVar<double> dLat; //NUMBER of cells in latitude
  GlobalVar<double> dLon; //NUMBER of cells in longitude
  GlobalVar<int> geodesic_l; //Geodesic grid level
  GlobalVar<double> period; // orbital period (=rotation period), assuming syncronous rotation.
  GlobalVar<double> endTime; // maximum simulation run time
  GlobalVar<std::string> potential; // string for tidal potential type
  GlobalVar<std::string> friction; // string for drag type
  GlobalVar<std::string> surface; // string for surface boundary condition
  GlobalVar<std::string> solver; // string for type of time sovler
  GlobalVar<bool> init; // boolean for using initial conditions
  GlobalVar<double> converge; // convergence criteria for ODIS.

  // Variables to switch on or off output of certain Variables
  // i.e., setting kinetic to true outputs global averaged kinetic energy
  GlobalVar<bool> diss; // dissipated energy
  GlobalVar<bool> kinetic; // kinetic energy
  GlobalVar<bool> work; // work flux - currently not implemented.

  GlobalVar<bool> field_displacement_output;
  GlobalVar<bool> field_velocity_output;
  GlobalVar<bool> field_diss_output;
  GlobalVar<bool> sh_coeff_output;

  // Time in fraction of orbital period for output of all parameters
  GlobalVar<double> outputTime;

  // constructor to initialise and/or read all variables from input file.
  Globals();
  Globals(int action); // action is either 1 (use defaults) or 0 (read input)

  // Member function to read global variables from input.in file.
  int ReadGlobals(void);
};



#endif
