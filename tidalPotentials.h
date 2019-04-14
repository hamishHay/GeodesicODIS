/* Contains functions to compute the components of the tidal potential gradient.
*
* Each function requires two Field objects to write the gradient in latitude
* and longitude, a Globals object containing constants relevant to the
* calculation, and the current simulation time.
*/

#ifndef TIDALPOTENTIALS_H
#define TIDALPOTENTIALS_H

#include "mesh.h"
#include "array2d.h"
#include <math.h>


// Degree 2 component of the eccentricity tide (see Tyler 2011, Matsuyama 2014)
void deg2Ecc(Mesh * grid, Array2D<double> &, double simulationTime, double radius, double omega, double ecc);

// Degree 2 component of the eccentricity libration tide (see Tyler 2011, Matsuyama 2014)
void deg2EccLib(Mesh * grid, Array2D<double> &, double simulationTime, double radius, double omega, double ecc);

// Degree 2 component of the eccentricity tide (see Tyler 2011, Matsuyama 2014)
void deg2EccWest(Mesh * grid, Array2D<double> &, double simulationTime, double radius, double omega, double ecc);

// Degree 2 component of the eccentricity tide (see Tyler 2011, Matsuyama 2014)
void deg2EccEast(Mesh * grid, Array2D<double> &, double simulationTime, double radius, double omega, double ecc);

// Degree 2 component of the eccentricity-radial tide (see Tyler 2011, Matsuyama 2014)
void deg2EccRad(Mesh * grid, Array2D<double> &, double simulationTime, double radius, double omega, double ecc);

// Degree 2 component of the obliquity tide (see Tyler 2011, Matsuyama 2014)
void deg2Obliq(Mesh * grid, Array2D<double> &, double simulationTime, double radius, double omega, double theta);

// Degree 2 component of the westward moving obliquity tide (see Tyler 2011, Matsuyama 2014)
void deg2ObliqWest(Mesh * grid, Array2D<double> &, double simulationTime, double radius, double omega, double theta);

// Degree 2 component of the eastward moving obliquity tide (see Tyler 2011, Matsuyama 2014)
void deg2ObliqEast(Mesh * grid, Array2D<double> &, double simulationTime, double radius, double omega, double theta);

void deg2Full(Mesh * grid, Array2D<double> &, double simulationTime, double radius, double omega, double theta, double ecc);

void deg2EccTotal(Mesh * grid, Array2D<double> &, double simulationTime, double radius, double omega, double theta, double ecc);

void deg2Planet(Mesh * grid, Array2D<double> &, Array1D<double> &, double simulationTime, double radius);

void deg2General(Mesh * grid, Array2D<double> &, double simulationTime, double radius, Array1D<double> & scalar);

void deg2PlanetObl(Mesh * grid, Array2D<double> &, double simulationTime, double radius, double obl);

//
// // Degree 2 component of the eccentricity-libration tide (see Tyler 2011, Matsuyama 2014)
// void deg2EccLib(Field * DUlat, Field * DUlon, double simulationTime, double radius, double omega, double ecc);
//
// // Degree 2 component of the obliquity tide (see Tyler 2011, Matsuyama 2014)
// void deg2Obliq(Field * dUlat, Field * dUlon, double simulationTime, double radius, double omega, double theta);
//
// // Degree 3 component of the eccentricity tide (see docs)
// void deg3Ecc(Field * DUlat, Field * DUlon, double simulationTime, double radius, double smAxis, double omega, double ecc);
//
// // Degree 3 component of the obliquity tide (see docs)
// void deg3Obliq(Field * DUlat, Field * DUlon, double simulationTime, double radius, double smAxis, double omega, double theta);


#endif
