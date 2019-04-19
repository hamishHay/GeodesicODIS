#ifndef INITIALCONDITIONS_H
#define INITIALCONDITIONS_H

#include "mesh.h"
#include "globals.h"
#include "array3d.h"
#include "array2d.h"
#include "array1d.h"

int loadInitialConditions(Globals *, Mesh * mesh,
                          Array2D<double> &, Array3D<double> &,
                          Array1D<double> &, Array2D<double> &);

int writeInitialConditions(Globals *, Mesh *,
                           Array2D<double> &, Array3D<double> &,
                           Array1D<double> &, Array2D<double> &);

#endif
