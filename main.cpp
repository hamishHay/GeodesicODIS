// SearsTidalDissipation.cpp : Defines the entry point for the console application.
//

#include "iostream"
#include "conio.h"
#include "mesh.h"
#include "field.h"
#include "globals.h"
#include "solver.h"

using namespace std;

void main(void)
{
	// Initialise and Read in Constants/Select options
	// These input parameters are constant in both space and time
	// and are thus declared as global variables
	//initialiseGlobals();

	// Initialise domain/read in initial condition
	//initialiseMesh();

	Globals * val = new Globals(1);
	//val->timeStep = 400;
	Mesh * grid = new Mesh(val); //Pass in a pointer to globals instance, and grid using dLat and dLon values.
	Field * u = new Field(grid,0,1); //Construct velocity storage field based around grid
	Field * v = new Field(grid,1,0); //Note velocity is staggered and lies on seperate nodes to eta and U
	Field * eta = new Field(grid,0,0);
	Field * dUlat = new Field(grid,1,0);
	Field * dUlon = new Field(grid,0,1);

	// Enter solver/time loop
	Solver * solution = new Solver(0, 1, val, grid, dUlon, dUlat,  u, v, eta);
	solution->Solve();
	// Simulation end

	_getch();

	// Unallocate all memory
	delete solution;
	delete val;
	delete grid;
	delete u;
	delete v;
	delete eta;
	delete dUlat;
	delete dUlon;
}

