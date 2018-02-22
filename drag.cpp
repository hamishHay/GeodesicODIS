#include "drag.h"
#include <math.h>
#include <iostream>
#include <sstream>


void linearDrag(int & node_num, double & drag_coeff, Array2D<double> & dv_dt, Array2D<double> & v)
{
    int i;

    for (i=0; i<node_num; ++i)
    {
        dv_dt(i,0) -= v(i,0) * drag_coeff;
        dv_dt(i,1) -= v(i,1) * drag_coeff;
    }
};


void quadraticDrag(int & node_num, double & drag_coeff, double & h, Array2D<double> & dv_dt, Array2D<double> & vel)
{
    int i;
    double inv_h;
    double u, v, vv, uu, sq_uv;

    inv_h = drag_coeff/h;

    for (i=0; i<node_num; ++i)
    {
      u = vel(i,0);
      v = vel(i,1);
      uu = u*u;
      vv = v*v;
      sq_uv = sqrt(uu + vv);

      dv_dt(i,0) -= inv_h * sq_uv * u;
      dv_dt(i,1) -= inv_h * sq_uv * v;
    }
};

void quadraticDrag2Layer(int & node_num, double & drag_coeff, double & den_ratio, double & radius_ratio, Array1D<double> & h_l1, Array1D<double> & h_l2, Array2D<double> & vel_shear, Array2D<double> & dv_dt)
{
    int i;
    double gam_x;
    double h_factor;
    double h1, h2;
    double u, v, vv, uu, sq_uv;

    gam_x = den_ratio*radius_ratio;

    for (i=0; i<node_num; i++)
    {
        h1 = h_l1(i);
        h2 = h_l2(i);

        h_factor = gam_x * (h1*h1*h1) + (h2*h2*h2) /
                   (h1*h2 * (gam_x*h1 + h2)*(gam_x*h1 + h2));

       u = vel_shear(i,0);
       v = vel_shear(i,1);
       uu = u*u;
       vv = v*v;
       sq_uv = sqrt(uu + vv);

       dv_dt(i, 0) -= drag_coeff * h_factor * sq_uv * u;
       dv_dt(i, 1) -= drag_coeff * h_factor * sq_uv * v;

    }

}
