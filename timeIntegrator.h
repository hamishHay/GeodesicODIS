#ifndef TIMEINTEGRATOR_H_INCDLUDED
#define TIMEINTEGRATOR_H_INCDLUDED

#include "mesh.h"
#include "globals.h"
#include <math.h>
#include <iostream>

int eulerIntegrator(Globals *, Mesh *);
int ab3Integrator(Globals *, Mesh *);

#endif
