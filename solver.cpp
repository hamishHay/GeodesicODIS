#include "mesh.h"
#include "globals.h"
#include "outFiles.h"
#include "solver.h"
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include <math.h>
#include <iostream>
#include <sstream>

std::ostringstream outstring;

int solveODIS(Globals * globals, Mesh * mesh)
{
    OutFiles * Output;
    Output = globals->Output;

    outstring << "Identifying using selected solver method...";

    

    Output->Write(OUT_MESSAGE, &outstring);

    return 1;
};
