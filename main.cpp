#include "mesh.h"
#include "field.h"
#include "globals.h"
#include "solver.h"
#include "mass.h"
#include "energy.h"

int main(void)
{
	//Read constants values (0) or use defaults (1)
	Globals * constants = new Globals(0);

	Mesh * grid = new Mesh(constants); //Pass in a pointer to globals instance, and grid using dLat and dLon values.
	Field * u = new Field(grid,0,1); //Construct velocity storage field based around grid
	Field * v = new Field(grid,1,0); //Note velocity is staggered and lies on seperate nodes to eta and U
	Field * eta = new Field(grid,0,0);
	Field * dUlat = new Field(grid,1,0);
	Field * dUlon = new Field(grid,0,1);

	//conserved quantities:
	Mass * mass = new Mass(grid, 0, 0, constants, eta);
	Energy * energy = new Energy(grid, 0, 0, constants, u, v, mass);

	// Enter solver/time loop
	Solver * solution = new Solver(0, 1, constants, grid, dUlon, dUlat,  u, v, eta, energy);
	solution->Solve();
	// Simulation end


	// Unallocate all memory
	delete solution;
	delete constants;
	delete grid;
	delete u;
	delete v;
	delete eta;
	delete dUlat;
	delete dUlon;

    return 0;
}
