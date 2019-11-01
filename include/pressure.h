#ifndef PRESSURE_H_INCDLUDED
#define PRESSURE_H_INCDLUDED

#include "mesh.h"
#include "globals.h"
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include "sphericalHarmonics.h"
#include <math.h>
#include <iostream>

int updatePressure(Globals *, Mesh *, Array1D<double> &, Array2D<double> &, Array2D<double> &);

#endif
