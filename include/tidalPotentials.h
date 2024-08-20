/* Contains functions to compute the components of the tidal potential gradient.
*
* Each function requires two Field objects to write the gradient in latitude
* and longitude, a Globals object containing constants relevant to the
* calculation, and the current simulation time.
*/

#ifndef TIDALPOTENTIALS_H
#define TIDALPOTENTIALS_H

#include "mesh.h"
#include "globals.h"
#include "array2d.h"
#include <math.h>

void forcing(Globals * consts, Mesh * grid, Array1D<double> & potential, int forcing_type, double time, double ecc=0.1, double obl=0.1);

// WOMA function 
// Get initial position and velocity of the impactor with initial distance r 
void impact_pos_vel_b_v_c_r(double &x, double &y, double &t_out_, double r, double b, double vC, double Rt, double Ri, double Mt, double Mi);

// WOMA function 
// Get the position of the impactor at time t
void impact_pos_vel_b_v_c_t(double &x, double &y, double simulationTime, double b, double vC, double Rt, double Ri, double Mt, double Mi, double rMaxFactor=100.0);

#endif

