#ifndef GLOBALS_H
#define GLOBALS_H

#include <string>
#include <array>
#include <vector>
#include "globalVar.h"
//#include "boost/variant.hpp"

extern double pi;
extern double radConv;

//Class stores all global values
class Globals {
private:
	void SetDefault(void);
public:
	std::vector<IGlobalVar *> allGlobals;
	
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
	GlobalVar<double> dLat;
	GlobalVar<double> dLon;
	GlobalVar<double> period;
	GlobalVar<int> endTime;
	GlobalVar<std::string> potential;

		
	Globals();
	Globals(int action);
	int ReadGlobals(void);
};



#endif