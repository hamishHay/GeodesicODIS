#include "energy.h"
#include "array1d.h"
#include "array2d.h"
#include "globals.h"
#include <math.h>
#include <iostream>
#include <sstream>


void updateEnergy(Globals * globals, double & avg_diss, Array1D<double> & e_diss, Array2D<double> & v, Array1D<double> & mass)
{
    int node_num, i;
    double drag_coeff, h, r;

    drag_coeff = globals->alpha.Value();
    h = globals->h.Value();
    r = globals->radius.Value();

    node_num = globals->node_num;

    avg_diss = 0.0;

    // UPDATE ENERGY
    switch (globals->fric_type)
    {
      case LINEAR:
          for (i = 0; i<node_num; i++)
          {
            e_diss(i) = drag_coeff * mass(i)
                            * (v(i,0)*v(i,0) + v(i,1)*v(i,1));

            avg_diss += e_diss(i);
        }
        break;

      case QUADRATIC:
          for (i = 0; i<node_num; i++)
          {
            e_diss(i) = drag_coeff/h * mass(i)
                                * sqrt(v(i,0)*v(i,0) + v(i,1)*v(i,1))
                                * (v(i,0)*v(i,0) + v(i,1)*v(i,1));
            avg_diss += e_diss(i);
        }
        break;
    }

    avg_diss /= (4. * pi * r * r);

};
