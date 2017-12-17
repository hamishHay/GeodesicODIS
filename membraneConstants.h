#ifndef MEMBRANECONSTANTS_H
#define MEMBRANECONSTANTS_H

#include "mesh.h"
#include "globals.h"
#include "array3d.h"
#include "array2d.h"
#include "array1d.h"

void membraneNuBeta(Array1D<double> &, Array1D<double> &, int l_max, Globals * consts);

void loveNumTide(Array1D<double> &, Array1D<double> &, int l_max, Globals * consts);
void loveNumLoad(Array1D<double> &, Array1D<double> &, int l_max, Globals * consts);

void effectiveRigidity(Array1D<double>, Globals * consts);

void springConst(Array1D<double> &, Array1D<double> &, int l_max, Globals * consts);

#endif
