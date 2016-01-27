#ifdef _WIN32
#include <direct.h>
#include <Windows.h>
#define getcwd _getcwd
#define SEP "\\"

#elif _WIN64
#include <direct.h>
#define getcwd _getcwd
#define SEP "\\"

#elif __linux__
#include <unistd.h>
#define SEP "/"

#else
#error "OS not supported!"
#endif

#include "globals.h"
#include "outFiles.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <math.h>

Globals::Globals() :Globals(1) {};

Globals::Globals(int action) {
	path = Output.path;

	strcpy(cpath, path.c_str());

	if (action) SetDefault();
	else {
		SetDefault();

		//Append variable IDs to all global variables, then add them to variable list
		//All variables described in globals.h
		radius.SetStringID("radius");
		allGlobals.push_back(&radius);

		angVel.SetStringID("angular velocity");
		allGlobals.push_back(&angVel);

		loveK2.SetStringID("k2");
		allGlobals.push_back(&loveK2);

		loveH2.SetStringID("h2");
		allGlobals.push_back(&loveH2);

		loveReduct.SetStringID("love reduction factor");
		allGlobals.push_back(&loveReduct);

		h.SetStringID("ocean thickness");
		allGlobals.push_back(&h);

		g.SetStringID("surface gravity");
		allGlobals.push_back(&g);

		a.SetStringID("semimajor axis");
		allGlobals.push_back(&a);

		e.SetStringID("eccentricity");
		allGlobals.push_back(&e);

		theta.SetStringID("obliquity");
		allGlobals.push_back(&theta);

		timeStep.SetStringID("time step");
		allGlobals.push_back(&timeStep);

		converge.SetStringID("converge");
		allGlobals.push_back(&converge);

		alpha.SetStringID("friction coefficient");
		allGlobals.push_back(&alpha);

		dLat.SetStringID("latitude spacing");
		allGlobals.push_back(&dLat);

		dLon.SetStringID("longitude spacing");
		allGlobals.push_back(&dLon);

		period.SetStringID("orbital period");
		allGlobals.push_back(&period);

		endTime.SetStringID("simulation end time");
		allGlobals.push_back(&endTime);

		potential.SetStringID("potential");
		allGlobals.push_back(&potential);

		friction.SetStringID("friction type");
		allGlobals.push_back(&friction);

		init.SetStringID("initial conditions");
		allGlobals.push_back(&init);

		diss.SetStringID("dissipation output");
		allGlobals.push_back(&diss);

		kinetic.SetStringID("kinetic output");
		allGlobals.push_back(&kinetic);

		work.SetStringID("work flux output");
		allGlobals.push_back(&work);

		outputTime.SetStringID("output time");
		allGlobals.push_back(&outputTime);

		ReadGlobals(); //Read globals from input.in file
	}

	OutputConsts();

	endTime.SetValue(endTime.Value()*period.Value());

	if (friction.Value() == "LINEAR") fric_type = LINEAR;
	else if (friction.Value() == "QUADRATIC") fric_type = QUADRATIC;
	else {
		outstring << "No friction type found." << std::endl;
		Output.Write(ERR_MESSAGE, &outstring);
		Output.TerminateODIS();
	}
};

int Globals::ReadGlobals(void) {
	//Function to read model input parameters from input.in only.
	//Avoids the need to recompile for each different model setup.
	std::ifstream inputFile(path + SEP + "input.in",std::ifstream::in);
	std::string line, varStr, varVal;

	int valInt;
	double valDouble;
	bool valBool = true;
	std::string valStr;
	bool added=false;

	if (inputFile.is_open())
	{

		outstring << "Found input file." << std::endl;
		Output.Write(OUT_MESSAGE, &outstring);

		while (std::getline(inputFile >> std::ws, line, ';')) { //>> std::ws removes any leading whitespace

			//Convert string to lower case
			std::transform(line.begin(), line.end(), line.begin(), ::tolower);
			std::getline(inputFile >> std::ws, varVal, ';');
			std::stringstream value(varVal);

			added = false;
			unsigned int i = 0; //unsigned to remove warnings
			do {
				if (line == allGlobals[i]->StringID()) {
					if (allGlobals[i]->IsType("double")) {
						value >> valDouble;
						((GlobalVar<double > *) allGlobals[i])->SetValue(valDouble);
						added = true;
						allGlobals[i]->Added(added);
					}
					else if (allGlobals[i]->IsType("int")) {
						value >> valInt;
						((GlobalVar<int> *) allGlobals[i])->SetValue(valInt);
						added = true;
						allGlobals[i]->Added(added);
					}
					else if (allGlobals[i]->IsType("bool")) {
						value >> valStr;
						if (valStr == "false") valBool = false;
						else if (varStr == "true") valBool = true;
						((GlobalVar<bool> *) allGlobals[i])->SetValue(valBool);
						added = true;
						allGlobals[i]->Added(added);
					}
					else if (allGlobals[i]->IsType("string")) {
						value >> valStr;
						((GlobalVar<std::string> *) allGlobals[i])->SetValue(valStr);
						added = true;
						allGlobals[i]->Added(added);
					}
					else {
						outstring << "Unsused type read for " << allGlobals[i]->StringID() << ". " << std::endl;
						Output.Write(ERR_MESSAGE, &outstring);
						Output.TerminateODIS();
					};
				}
				i++;
			} while (i < allGlobals.size() && !added);

			std::getline(inputFile, line, ';'); //Ignores last part of input file line
		};

		inputFile.close();



		for (unsigned int i = 0; i < allGlobals.size(); i++){
			if (!allGlobals[i]->IsAdded()) {
				outstring << "WARNING: Unassigned global constant: " << allGlobals[i]->StringID() << ". " << std::endl;
				outstring << "Assigning default value could be risky..." << std::endl;
				outstring << "Using Titan value of ";

				if (allGlobals[i]->IsType("double")) {
					outstring << ((GlobalVar<double > *) allGlobals[i])->Value() << std::endl << std::endl;
				}
				else if (allGlobals[i]->IsType("int")) {
					outstring << ((GlobalVar<int > *) allGlobals[i])->Value() << std::endl << std::endl;
				}
				else if (allGlobals[i]->IsType("bool")) {
					outstring << ((GlobalVar<bool > *) allGlobals[i])->Value() << std::endl << std::endl;
				}
				else if (allGlobals[i]->IsType("string")) {
					outstring << ((GlobalVar<std::string > *) allGlobals[i])->Value() << std::endl << std::endl;
				}

				Output.Write(ERR_MESSAGE, &outstring);
			}
		}
	}
	else {

		outstring << "Unable to open 'input.in' file." << std::endl << std::endl;
		Output.Write(ERR_MESSAGE, &outstring);
		Output.TerminateODIS();
	}
	return 0;
};

void Globals::SetDefault(void) {

	//Satellite Radius
	radius.SetValue(2574.73e3);

	//k_2 love number
	loveK2.SetValue(.120);

	//h_2 love number
	loveH2.SetValue(.2);

	//Love reduction factor
	//loveReduct.SetValue(1 + loveK2.Value() - loveH2.Value());
	loveReduct.SetValue(0.920012);
	//Ocean thickness
	h.SetValue(400);

	//Surface gravity
	g.SetValue(1.35);

	//Semimajor axis
	a.SetValue(1221.87e6);

	//Eccentricity
	e.SetValue(0.0288);

	//Obliquity
	theta.SetValue(0.32*pi / 180.);

	//Coefficient of linear friction
	alpha.SetValue(2.28e-7);

	//Lat and long spacing
	dLat.SetValue(45);
	dLon.SetValue(90);

	//Orbital period
	period.SetValue(16. * 24 * 60 * 60);

	//Time step
	timeStep.SetValue(40.);

	//Convergence Criteria
	converge.SetValue(1e-7);

	//Angular velocity
	angVel.SetValue(4.56e-6);

	//Orbital period
	period.SetValue(2*pi/angVel.Value());

	//Simulation end time
	endTime.SetValue(80.);

	//Potential
	potential.SetValue("ECC");

	//Friction Type
	friction.SetValue("QUADRATIC");

	//Initial conditions
	init.SetValue(false);

	//Output dissipated energy?
	diss.SetValue(true);

	//Output kinetic energy?
	kinetic.SetValue(false);

	//Output work flux?
	work.SetValue(false);

	//Output time in fraction of orbital period. Default is set to output
	//every period.
	outputTime.SetValue(1);
};

void Globals::OutputConsts(void){
	// Print all constants to output file.
	outstring << "Model parameters: " << std::endl;
	for (unsigned int i = 0; i < allGlobals.size(); i++){

		outstring << allGlobals[i]->StringID() << ":\t\t\t\t ";
		if (allGlobals[i]->IsType("double")) {
			outstring << ((GlobalVar<double > *) allGlobals[i])->Value();
		}
		else if (allGlobals[i]->IsType("int")) {
			outstring << ((GlobalVar<int > *) allGlobals[i])->Value();
		}
		else if (allGlobals[i]->IsType("bool")) {
			outstring << ((GlobalVar<bool > *) allGlobals[i])->Value();
		}
		else if (allGlobals[i]->IsType("string")) {
			outstring << ((GlobalVar<std::string > *) allGlobals[i])->Value();
		}

		Output.Write(OUT_MESSAGE, &outstring);
	}
};
