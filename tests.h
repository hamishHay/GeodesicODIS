#ifndef TESTS_H_INCDLUDED
#define TESTS_H_INCDLUDED

#include "mesh.h"
#include "globals.h"
#include "array2d.h"
#include "array1d.h"
#include <math.h>
#include <iostream>

void runOperatorTests(Globals *, Mesh *);


void setBeta(Array1D<double> &, int, int, int, Array2D<double> &);

void setU(Array2D<double> &, int, int, int, Array2D<double> &);

void setDivU(Array1D<double> &, int, int, int, double, Array2D<double> &);

void setGradBeta(Array2D<double> &, int, int, int, double, Array2D<double> &);

#endif
