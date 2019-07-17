#include "energy.h"
#include "array1d.h"
#include "array2d.h"
#include "globals.h"
#include "interpolation.h"
#include <math.h>
#include <iostream>
#include <sstream>


void updateEnergy(Globals * globals, double & avg_flux, Array1D<double> & e_flux, Array2D<double> & vel)
{
    int node_num, face_num, i;
    double drag_coeff, h, r;
    // double avg_flux;

    drag_coeff = globals->alpha.Value();
    h = globals->h.Value();
    r = globals->radius.Value();

    node_num = globals->node_num;
    face_num = globals->face_num;

    avg_flux = 0.0;

    // UPDATE ENERGY
    switch (globals->fric_type)
    {
      case LINEAR:
          for (i = 0; i<face_num; i++)
          {
            e_flux(i) = drag_coeff * 1000.0 * h * (vel(i,0)*vel(i,0) + vel(i,1)*vel(i,1));

            avg_flux += e_flux(i);

            // e_flux(i) *= 1./mass(i) * 1000.0 * h; // Convert to flux?



            // flux += 1000.0 * h * drag_coeff * (v(i,0)*v(i,0) + v(i,1)*v(i,1));
        }
        break;

      case QUADRATIC:
          for (i = 0; i<face_num; i++)
          {
            e_flux(i) = drag_coeff/h * sqrt(vel(i,0)*vel(i,0) + vel(i,1)*vel(i,1))
                                * (vel(i,0)*vel(i,0) + vel(i,1)*vel(i,1));


            avg_flux += e_flux(i);
        }
        break;
    }

    // avg_flux /= (4. * pi * r * r);
    // TODO - convert to area weighted average
    avg_flux = avg_flux/face_num;

};
