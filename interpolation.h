#ifndef INTERPOLATION_H_INCDLUDED
#define INTERPOLATION_H_INCDLUDED

#include "mesh.h"
#include "globals.h"
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include <math.h>
#include <iostream>

int interpolateGG2LL(Globals *, Mesh *, Array2D<double> &, Array1D<double> &, Array3D<double> &, Array3D<double> &, Array2D<int> &);

#endif
