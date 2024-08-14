#ifndef MOMADVECTION_H
#define MOMADVECTION_H

#include "mesh.h"
#include "globals.h"
#include "array1d.h"

void calculateMomentumAdvection(Globals *globals,
                                Mesh *mesh,
                                Array1D<double> &dvdt,            // array to apply advection to
                                Array1D<double> &vel,             // face normal velocities
                                Array1D<double> &layer_thickness, // node ocean thickness
                                Array1D<double> &Ekin);

#endif