// FILE: viscosity.h
// DESCRIPTION: header file for the viscosity function. Viscosity calculates
//              the viscosity term for the two velocity fields u and v using
//              second order accurate central differences or first order forward
//              and backward differences around any boundaries. Polar values can
//              either be interpolated or calulated with foward/backward
//              differences.
//
//   Version              Date                   Programmer
//     0.1      |        19/05/2017        |        H. Hay        |
// ---- Adding viscosity term to ODIS for stability


#ifndef VISCOSITY_H
#define VISCOSITY_H

#include "globals.h"
#include "field.h"

#include <math.h>

void UpdateViscosity(Field * u, Field * v, Globals * constants);

#endif
