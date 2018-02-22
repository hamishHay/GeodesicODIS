#ifndef ENERGY_H_INCDLUDED
#define ENERGY_H_INCDLUDED

#include "mesh.h"
#include "globals.h"
#include "array2d.h"
#include <math.h>
#include <iostream>

void updateEnergy(Globals *, double &, Array1D<double> &, Array2D<double> &, Array1D<double> &);
void updateEnergy2Layer(Globals *, double &, Array1D<double> &, Array1D<double> &, Array2D<double> &, Array1D<double> &);
#endif
