#ifndef SPHERICALHARMONICS_H
#define SPHERICALHARMONICS_H

#include "mesh.h"
#include "globals.h"
#include "array3d.h"
#include "array2d.h"
#include "array1d.h"

void getSHCoeffs(Mesh *, Array1D<double> &, Array2D<double> &);

#endif
