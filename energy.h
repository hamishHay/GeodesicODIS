#ifndef ENERGY_H_INCDLUDED
#define ENERGY_H_INCDLUDED

#include "mesh.h"
#include "globals.h"
#include "array2d.h"
#include <math.h>
#include <iostream>

void updateEnergy(Globals *, double &, Array1D<double> &, Array2D<double> &, Array1D<double> &);

void calculateKE(Array1D<double> & ke, Array1D<double> & velocity, Array1D<double> & face_areas, Array1D<double> & node_areas, int face_num, int node_num);

#endif
