#ifndef GLOBALS_H
#define GLOBALS_H

#include <string>
#include <array>
#include <vector>
#include "globalVar.h"
#include "outFiles.h"
#include <Windows.h>

extern double pi;
extern double radConv;

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
	char cpath[MAX_PATH];
	
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
	GlobalVar<std::string> friction;
	GlobalVar<bool> init;
		
	Globals();
	Globals(int action);
	int ReadGlobals(void);
};



#endif