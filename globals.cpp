#include "globals.h"
#include "outFiles.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

double pi = 3.1415926535897932384626433832795028841971693993751058; //change for derrick
double radConv = pi / 180;

Globals::Globals() :Globals(1) {};

Globals::Globals(int action) {
	char buffer[MAX_PATH];
	GetModuleFileName(NULL, buffer, MAX_PATH);
	std::string::size_type pos = std::string(buffer).find_last_of("\\/");
	path = std::string(buffer).substr(0, pos);

	strcpy(cpath, path.c_str());

	if (action) SetDefault();
	else {
		SetDefault();

		//Append variable IDs;
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

		init.SetStringID("initial conditions");
		allGlobals.push_back(&init);

		ReadGlobals(); //Read globals from input.in file
	}

	OutputConsts();

	endTime.SetValue(endTime.Value()*period.Value());
};

int Globals::ReadGlobals(void) {
	//Function to read model input parameters from input.in only.
	//Avoids the need to recompile for each different model setup.
	std::ifstream inputFile(path + "\\input.in",std::ifstream::in);
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
			int i = 0;
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

			std::getline(inputFile, line, ';'); //Ignores last part of input file
		};

		inputFile.close();

		for (int i = 0; i < allGlobals.size(); i++){
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
				//TerminateODIS();
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
	radius.SetValue(2574.7e3);

	//k_2 love number
	loveK2.SetValue(.120);

	//h_2 love number
	loveH2.SetValue(.2);

	//Love reduction factor
	loveReduct.SetValue(1 + loveK2.Value() - loveH2.Value());

	//Ocean thickness
	h.SetValue(400);

	//Surface gravity
	g.SetValue(1.3);

	//Semimajor axis
	a.SetValue(1221865e3);

	//Eccentricity
	e.SetValue(0.0288);

	//Obliquity 
	theta.SetValue(0.32*pi / 180);

	//Coefficient of linear friction
	alpha.SetValue(2.28e-7);

	//Lat and long spacing
	dLat.SetValue(2);
	dLon.SetValue(2);

	//Orbital period
	period.SetValue(16 * 24 * 60 * 60);

	//Time step
	timeStep.SetValue(40.);

	//Convergence Criteria
	converge.SetValue(1e-7);

	//Angular velocity
	angVel.SetValue(2 * pi / period.Value());

	//Simulation end time
	endTime.SetValue(80);

	//Potential
	potential.SetValue("OBLIQ");

	//Initial conditions
	init.SetValue(false);
};

void Globals::OutputConsts(void){
	outstring << "Model parameters: " << std::endl;
	for (int i = 0; i < allGlobals.size(); i++){

		outstring << allGlobals[i]->StringID() << ": ";
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

		Output.Write(OUT_MESSAGE, &outstring);
	}
};
