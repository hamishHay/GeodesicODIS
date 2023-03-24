#ifndef UPDATEETA_H_INCDLUDED
#define UPDATEETA_H_INCDLUDED

#include "mesh.h"
#include "globals.h"
#include "array1d.h"

int updateEta(Globals *, 
              Mesh *, 
              Array1D<double> &detadt, 
              Array1D<double> &v_t0, 
              Array1D<double> &eta,
              Array1D<double> &h_total);

#endif