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

const double pi = 3.1415926535897932384626433832795028841971693993751058;
const double radConv = pi / 180;
const int PATH = 512;

enum Friction { LINEAR, QUADRATIC };

//Class stores all global values
class Globals {
private:
	void SetDefault(void);
	void OutputConsts(void);

public:
	OutFiles Output;
	std::vector<IGlobalVar *> allGlobals;

	Friction fric_type;

	std::ostringstream outstring;

	std::string path;
	char cpath[512];

	GlobalVar<double> angVel;
	GlobalVar<double> radius; //Object Radius
	GlobalVar<double> loveK2; // k_2 Love Number
	GlobalVar<double> loveH2; // h_2 Love Number
	GlobalVar<double> loveReduct; //Love's reduction factor
	GlobalVar<double> h; //Ocean thickness
	GlobalVar<double> g; //Surface Gravity
	GlobalVar<double> a; //SemiMajor Axis
	GlobalVar<double> e; //Eccentricity
	GlobalVar<double> theta; //obliquity
	GlobalVar<double> timeStep;
	GlobalVar<double> alpha;
	GlobarVar<int> l_max; //Maximum spherical harmonic degree for expansions
	GlobalVar<double> dLat; //NUMBER of cells in latitude
	GlobalVar<double> dLon; //NUMBER of cells in longitude
	GlobalVar<double> period;
	GlobalVar<double> endTime;
	GlobalVar<std::string> potential;
	GlobalVar<std::string> friction;
	GlobalVar<bool> init;
	GlobalVar<double> converge;

	// Variables to switch on or off output of certain Variables
	// i.e., setting kinetic to true outputs global averaged kinetic energy
	GlobalVar<bool> diss;
	GlobalVar<bool> kinetic;
	GlobalVar<bool> work;

	// Time in fraction of orbital period for output of all parameters
	GlobalVar<double> outputTime;

	Globals();
	Globals(int action);

	int ReadGlobals(void);
};



#endif
