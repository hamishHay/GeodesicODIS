#ifndef INTERPOLATION_H_INCDLUDED
#define INTERPOLATION_H_INCDLUDED

#include "mesh.h"
#include "globals.h"
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include <math.h>
#include <iostream>

int interpolateGG2LLConservative(Globals *, Mesh *, Array2D<double> &, Array1D<double> &);

int interpolateVelocity(Globals *, Mesh *, Array2D<double> &, Array1D<double> &);

int interpolateLSQFlux(Globals *, Mesh *, Array1D<double> &, Array1D<double> &, Array1D<double> &, double beta=0.5);

#endif
