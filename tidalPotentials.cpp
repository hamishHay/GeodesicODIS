#include "field.h"
#include "globals.h"
#include "tidalPotentials.h"
#include <math.h>

/* Degree 2 component of the obliquity tide (see Tyler 2011, Matsuyama 2014)
 *
 * Finds the gradient components of the degree-2 obliquity tide for Fields
 * DUlat and DUlon. Only current simulation required is absolute simulation
 * time. The function directly writes the solution arrays in both Field
 * classes.
 *
 *      inputs: Field * DUlat           Field for latitude gradient of potential 
 *              Field * DUlon           Field for longitude gradient of potential
 *              double simulationTime   Current time in the simulation
 *              double radius           Satellite radius
 *              double omega            Satellite rotational angular speed
 *              double theta            Satellite obliquity in radians
 *
*/
void deg2Obliquity(Field * DUlat, Field * DUlon, double simulationTime, double radius, double omega, double theta)
{       
    double cosM, factor;                    // cos(Mean anomaly)
    double * cosLat, * cos2Lat, * sinLat;   // Pointer to 1D arrays for trig functions of latitude
    double * cosLon, * sinLon;              // Pointer to 1D arrays for trug functions of longitude
    double ** latGrad, ** lonGrad;          // Pointer to 2D array solutions
    int i,j;
    int latLen, lonLen;

    factor = -3.*pow(omega,2.0)*pow(radius,2.0)*theta;
    cosM = cos(omega*simulationTime);

    // Get pointer to arrays
    latGrad = DUlat->solution;
    lonGrad = DUlon->solution;

    // Set variables for dUlat solution
    latLen = DUlat->fieldLatLen;
    lonLen = DUlat->fieldLonLen;

    cos2Lat = DUlat->cos2Lat;
    cosLon = DUlat->cosLon;

    // Solve for dUdlat
    for (i=0; i<latLen; i++) {
        for (j=0; j<lonLen; j++) {
            latGrad[i][j] = factor*cos2Lat[i]*cosLon[j]*cosM;
        }
    }

    // Set variable for dUlon solution
    latLen = DUlon->fieldLatLen;
    lonLen = DUlon->fieldLonLen;

    cosLat = DUlon->cosLat;
    sinLat = DUlon->sinLat;
    sinLon = DUlon->sinLon;

    // Solve for dUdlon
    for (i=0; i<latLen; i++) {
        for (j=0; j<lonLen; j++) {
            lonGrad[i][j] = -factor*cosLat[i]*sinLat[i]*sinLon[j]*cosM;
        }
    }

};
