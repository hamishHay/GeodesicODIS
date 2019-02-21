#include "mesh.h"
#include "globals.h"
#include "outFiles.h"
#include "solver.h"
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include "timeIntegrator.h"
#include <math.h>
#include <iostream>
#include <sstream>


// Function to identify time solver method and call the corresponding function.
int solveODIS(Globals * globals, Mesh * mesh)
{
    std::ostringstream outstring;
    OutFiles * Output;
    Output = globals->Output;

    outstring << "Identifying using selected solver method... ";

    int err;

    switch (globals->solver_type)
    {
    case EULER:
        outstring << globals->solver.Value() << std::endl << std::endl;
        outstring << "Entering time solver " << globals->solver.Value() << "...";
        outstring << std::endl << std::endl;
        Output->Write(OUT_MESSAGE, &outstring);

        eulerExplicit(globals, mesh);
        break;

    case AB3:
        outstring << globals->solver.Value() << std::endl << std::endl;
        outstring << "Entering time solver " << globals->solver.Value() << "...";
        outstring << std::endl << std::endl;
        Output->Write(OUT_MESSAGE, &outstring);

        err = ab3Explicit(globals, mesh);
        break;

    default:
        outstring << "ERROR: NO SOLVER FOUND!" <<std::endl << std::endl;
        Output->Write(ERR_MESSAGE, &outstring);
        Output->TerminateODIS();
        break;
    }

    if (err) {
        outstring << "SOLVER RETURNED WITH AN ERROR..." <<std::endl << std::endl;
        Output->Write(OUT_MESSAGE, &outstring);
    }
    else {
        outstring << "Calculations appear to have finished!" << std::endl;
        Output->Write(OUT_MESSAGE, &outstring);
    }
    return 1;
};
