// FILE: interpolation.h
// DESCRIPTION: header file for interpolating functions. These are used to
//              interpolate velocity and displacement fields onto neighbouring
//              nodes.
//
//   Version              Date                   Programmer
//     0.1      |        19/05/2017        |        H. Hay        |
// ---- Adding inpolation functions to aid advection calculations


#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "field.h"

#include <math.h>

void InterpolateVel2Vel(Field * u, Field * v);
void InterpolateVel2Disp(Field * u, Field * v, Field * eta);
void InterpolateDisp2Vel(Field * u, Field *v, Field * eta);

#endif
