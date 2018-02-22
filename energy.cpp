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
    double flux;

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
            // flux += 1000.0 * h * drag_coeff * (v(i,0)*v(i,0) + v(i,1)*v(i,1));
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
    // avg_diss = flux/node_num;

};

void updateEnergy2Layer(Globals * globals, double & avg_diss, Array1D<double> & e_diss, Array1D<double> & h_l1, Array2D<double> & vel_shear, Array1D<double> & mass)
{
    int node_num, i;
    double drag_coeff, den_ratio, radius_ratio;
    double gam_x;
    double H, h1, h2, h_factor, r;

    drag_coeff = globals->alpha.Value();
    H = globals->h.Value();
    r = globals->radius_top.Value();

    den_ratio = globals->den_ratio.Value();
    radius_ratio = globals->radius_ratio.Value();

    gam_x = den_ratio * radius_ratio;

    node_num = globals->node_num;

    avg_diss = 0.0;

    // UPDATE ENERGY
    switch (globals->fric_type)
    {
      case LINEAR:
        break;

      case QUADRATIC:
          for (i = 0; i<node_num; i++)
          {
              h1 = h_l1(i);
              h2 = H - h1;

              h_factor = (drag_coeff/H) * ((h1*h1*h1) + (h2*h2*h2))/((gam_x*h1 + h2)*(gam_x*h1 + h2)*(gam_x*h1 + h2));

            e_diss(i) = h_factor * mass(i)
                                * sqrt(vel_shear(i,0)*vel_shear(i,0) + vel_shear(i,1)*vel_shear(i,1))
                                * (vel_shear(i,0)*vel_shear(i,0) + vel_shear(i,1)*vel_shear(i,1));
            avg_diss += e_diss(i);
        }
        break;
    }

    avg_diss /= (4. * pi * r * r);

};
