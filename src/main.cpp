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
//#include "tidalPotentials.h"
#include "array2d.h"
#include "solver.h"
#include "tests.h"
#include "gridConstants.h"


#include <iostream>

int main(void)
{
  // Read in global constants from input.in file with 0, and use defaults
  // with 1 (currently set to Titan parameters). All constants are stored in the
  // class "Globals".

  Globals * constants = new Globals(0);

  // Create the numerical grid using the minimum node spacing from "constants"
  // Globals instance.

  Mesh * grid = new Mesh(constants, constants->node_num, constants->face_num, constants->vertex_num, (int)constants->dLat.Value(), constants->l_max.Value());

  constants->OutputConsts();

  #if defined(TEST_OPERATORS)
    runOperatorTests(constants, grid);
  #else
    solveODIS(constants, grid);
  #endif 

  return 0;
}
