// FILE: advection.h
// DESCRIPTION: header file for the advection function. Advection calculates the
//              advection term found in the shallow water equations. The
//              quantities that are advected are velocity and mass (displacement).
//              Advection is performed using the upwind finite difference scheme.
//
//   Version              Date                   Programmer
//     0.1      |        24/05/2017        |        H. Hay        |
// ---- Adding advection term for high velocity simulations


#ifndef ADVECTION_H
#define ADVECTION_H

#include "globals.h"
#include "field.h"

#include <math.h>

void UpdateAdvection(Field * u, Field * v, Field * eta, Globals * constants);

#endif
