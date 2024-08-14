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

int interpolateVelocityCartRBF(Globals * globals, Mesh * mesh, Array2D<double> & interp_vel_xyz, Array1D<double> & normal_vel);

int interpolateLSQFlux(Globals *, Mesh *, Array1D<double> &, Array1D<double> &, Array1D<double> &, double beta=1.0);

#endif
