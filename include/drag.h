#ifndef DRAG_H_INCDLUDED
#define DRAG_H_INCDLUDED

#include "mesh.h"
#include "globals.h"
#include "array2d.h"
#include <math.h>
#include <iostream>

void linearDrag(int &, double &, Array2D<double> &, Array2D<double> &);
void quadraticDrag(int &, double, double, Array1D<double> &, Array2D<double> &, Array1D<double> &);

#endif
