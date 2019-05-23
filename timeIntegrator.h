#ifndef TIMEINTEGRATOR_H_INCDLUDED
#define TIMEINTEGRATOR_H_INCDLUDED

#include "mesh.h"
#include "globals.h"
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include <math.h>
#include <iostream>

int ab3Explicit(Globals *, Mesh *);
int rk4Explicit(Globals *, Mesh *);

int updateVelocity(Globals *, Mesh *, Array2D<double> &, Array2D<double> &, Array1D<double> &, double, int);
int updateDisplacement(Globals *, Mesh *, Array1D<double> &, Array2D<double> &, Array1D<double> &);

#endif
