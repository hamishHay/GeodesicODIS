#ifndef SPHERICALHARMONICS_H
#define SPHERICALHARMONICS_H

#include "mesh.h"
#include "globals.h"
#include "array3d.h"
#include "array2d.h"
#include "array1d.h"

void getSHCoeffsGG(Array2D<double> &, Array1D<double> &, Array3D<double> &, int, int);
void getSHCoeffsLL(Array2D<double> &, Array3D<double> &, int, int);
void getLegendreFuncs(double, Array2D<double> &, int);
void getLegendreFuncsDeriv(double, Array2D<double> &, int);

#endif
