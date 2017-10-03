#ifndef SPATIALOPERATORS_H
#define SPATIALOPERATORS_H

#include "mesh.h"
#include "globals.h"
#include "array3d.h"
#include "array2d.h"
#include "array1d.h"

void pressureGradient(Mesh *, Array2D<double> &, Array1D<double> &, double);

void pressureGradientSH(Globals *, Mesh *, Array2D<double> &, Array1D<double> &, double);

void velocityDivergence(Mesh *, Array1D<double> &, Array2D<double> &, double &, double);

void velocityDiffusion(Mesh *, Array2D<double> &, Array2D<double> &, double);

void scalarDiffusion(Mesh *, Array1D<double> &, Array1D<double> &, double);

#endif
