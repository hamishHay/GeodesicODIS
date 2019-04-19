#ifndef TEMPORALOPERATORS_H_INCDLUDED
#define TEMPORALOPERATORS_H_INCDLUDED

#include "mesh.h"
#include "globals.h"
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include <math.h>
#include <iostream>

int integrateAB3vector(Globals *, Mesh *, Array2D<double> &, Array3D<double> &, int interation);

int integrateAB3scalar(Globals *, Mesh *, Array1D<double> &, Array2D<double> &, int interation);

#endif
