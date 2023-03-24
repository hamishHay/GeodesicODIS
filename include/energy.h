#ifndef ENERGY_H_INCDLUDED
#define ENERGY_H_INCDLUDED

#include "mesh.h"
#include "globals.h"
#include "array2d.h"
#include <math.h>
#include <iostream>

void updateEnergy(Globals *, double &, Array1D<double> &, Array2D<double> &, Array1D<double> &);

void calculateKE(Globals *, Mesh *, Array1D<double> & ke, Array1D<double> & velocity);

#endif
