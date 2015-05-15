#ifndef GLOBALS_H
#define GLOBALS_H

#include <string>
#include <array>
#include <vector>
#include "globalVar.h"
//#include "boost/variant.hpp"

extern double pi;
/*
template<class TYPE>
class GlobalVar {
public:

	std::string stringID;
	TYPE inputValue;
	bool isEnd = false;

	GlobalVar * first;
	GlobalVar * next;

	GlobalVar();
	GlobalVar(std::string);

	TYPE Val(void);
};*/

//Class stores all global values
class Globals {
public:
	/*
	double angVel;
	double radius; //Object Radius
	double loveK2; // k_2 Love Number
	double loveH2; // h_2 Love Number
	double loveReduct; //Love's reduction factor
	double h; //Ocean thickness
	double g; //Surface Gravity
	double a; //SemiMajor Axis
	double e; //Eccentricity
	double theta; //obliquity
	double timeStep;
	double alpha;
	double dLat;
	double dLon;
	double period;
	int endTime; //End time, in unit of orbits
	*/

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

		
	Globals();
	Globals(int action);
	int ReadGlobals(void);
};



#endif