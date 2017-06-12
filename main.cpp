/* Welcome to Ocean Dissipation in Icy Satellites, a finite difference
 * computational fluid dynamics code for simulating global scale ocean dynamics
 * for any planetary body with a global ocean.
 *
 * The intended purpose of this code is to provide a high performance method of
 * quantifying ocean tidal dissipation, although the code can be modified to
 * handle any forcing.
 *
 * ODIS is also a useful tool for understanding the ocean dynamics of planetary
 * bodies, including in icy satellite and magma oceans. It is my hope that
 * the planetary community may find this open source code useful, but it should
 * also be pointed out that it is in development by only myself, someone who has
 * little rigorous training in numerical methods. As such, this version of ODIS
 * is intentionally constructed using basic numerical approaches. These include
 * the approximation of spatial derivatives using a finite difference method,
 * temporal derivatives using a first order Euler scheme, and a fixed grid in
 * latitde and longitude.
 *
 * FILE: main.cpp
 * DESCRIPTION: main program function. Controls allocation of all solvable and
 *              and transportable quantities in the model domain, as well as
 *              global constants and generation of the numerical grid.
 *
 *   Version              Date                 Programmer
 *      0.1     |       29/07/2016       |        H. Hay        |
 * ---- Initial version of ODIS, described in Hay and Matsuyama (2017)
*/

#include "outFiles.h"
#include "mesh.h"
// #include "field.h"
// #include "depth.h"
#include "globals.h"
// #include "solver.h"
// #include "mass.h"
// #include "energy.h"
#include "tidalPotentials.h"
#include "array2d.h"
#include "solver.h"


#include <iostream>

int main(void)
{
  // Read in global constants from input.in file with 0, and use defaults
  // with 1 (currently set to Titan parameters). All constants are stored in the
  // class "Globals".

  Globals * constants = new Globals(0);

  // Create the numerical grid using the minimum node spacing from "constants"
  // Globals instance.

  Mesh * grid = new Mesh(constants, constants->node_num);

  solveODIS(constants, grid);

  // Array2D<double> * velocity = new Array2D<double>(constants->node_num, 2);
  //
  //
  // double ** p;
  // p = new double *[constants->out_tags.size()-1];
  //
  // double time = 0.0;
  // double end_time = constants->endTime.Value();
  // double dt = 0.1 * end_time;
  // int count = 1;
  //
  //
  // while (time <= end_time)
  // {
  //     deg2Ecc(grid, velocity, time, 252.1e3, 5.31e-5, 0.0047);
  //
  //     p[0] = &(*velocity)(0,0);
  //
  //     constants->Output->DumpData(constants, count, p);
  //
  //     count++;
  //     time += dt;
  //
  // }


  // Allocate the solvable "Fields". These are all the quantities which ODIS aims
  // to calculate a solution for (velocity and surface displacement), as well as
  // the tidal forcing components. Each field is allocated to a particular node
  // in the grid, and is either staggered one node East or one node South of the
  // parent cell.

  // Field * u = new Field(grid,0,1); // Eastward velocity component. Staggered East.
  // Field * v = new Field(grid,1,0); // Northward velocity component. Staggered south
  // Field * eta = new Field(grid,0,0); // Surface displacement. Cell centered.
  // Field * dUlat = new Field(grid,1,0); // Latitudinal tidal potential gradient. Staggered south.
  // Field * dUlon = new Field(grid,0,1); // Longitudinal tidal potential gradient. Staggered east.

  // Depth * h = new Depth(grid);

  // Allocate memory for (hopefully) conserved quantities, mass and energy.
  // Classes Mass and Energy use Field as a parent class, and must be passed the
  // quantities on which they are derived.

  // Mass * mass = new Mass(grid, 0, 0, constants, eta, h);
  // Energy * energy = new Energy(grid, 0, 0, constants, u, v, mass);

  // Allocate space for class tupe Solver. Solver creates the environment to
  // begin the numerical calculations, with pointer access to all necessary
  // quantiies.

  // Solver * solution = new Solver(0, 1, constants, grid, dUlon, dUlat,  u, v, eta, energy, h);

  // Begin the calculations by calling the Solve member function. Probably
  // overkill.
  // solution->Solve();

  // Calculations are finished or ODIS has terminated with an error. Return from
  // main.

  return 0;
}
