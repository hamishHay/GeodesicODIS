#ifndef UPDATEMOMENTUM_H_INCDLUDED
#define UPDATEMOMENTUM_H_INCDLUDED

#include "array1d.h"
#include "globals.h"
#include "mesh.h"

int updateMomentum(Globals * constants,
                   Mesh * grid, 
                   Array1D<double> & dvdt, 
                   Array1D<double> & v_tm1, 
                   Array1D<double> & p_tm1, 
                   Array1D<double> & h_total, 
                   Array1D<double> & ekin, 
                   double GAMMA=0.0, 
                   double IMPLICIT=0.0);

#endif