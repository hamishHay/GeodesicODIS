#ifndef ANALYTICALLTE_H
#define ANALYTICALLTE_H

#include "mesh.h"
#include "globals.h"


std::array<double, 6> analyticalLTE(Globals * consts, Mesh * grid, int forcing_type, double lat, double lon, double t);

#endif