#ifndef TIDALPOTENTIALS_H
#define TIDALPOTENTIALS_H

#include "field.h"
#include "globals.h"
#include <math.h>

/* Contains functions to compute the components of the tidal potential gradient.
 *
 * Each function requires two Field objects to write the gradient in latitude
 * and longitude, a Globals object containing constants relevant to the
 * calculation, and the current simulation time.
*/

// Degree 2 component of the obliquity tide (see Tyler 2011, Matsuyama 2014)
void deg2Obliquity(Field * dUlat, Field * dUlon, double simulationTime, double radius, double angVel, double theta);

#endif
