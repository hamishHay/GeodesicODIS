#include "field.h"
#include "globals.h"
#include "tidalPotentials.h"
#include <math.h>

/* Degree 2 component of the eccentricity tide (see docs)
 *
 * Finds the gradient components of the degree-3 eccentricity tide for Fields
 * DUlat and DUlon. Only current simulation required is absolute simulation
 * time. The function directly writes the solution arrays in both Field
 * classes. The eccentricity tide contains terms for both P_20 and P_22
 * associated Legendre polynomials.
 *
 *      inputs: Field * DUlat           Field for latitude gradient of potential
 *              Field * DUlon           Field for longitude gradient of potential
 *              double simulationTime   Current time in the simulation
 *              double radius           Satellite radius
 *              double omega            Satellite rotational angular speed
 *              double ecc              Satellite eccentricity (must be << 1)
 *
*/
void deg2Ecc(Field * DUlat, Field * DUlon, double simulationTime, double radius, double omega, double ecc) {
    double cosM, sinM, factor;              // cos(Mean anomaly), sin(Mean anomaly)
    double * sin2Lat, * cosSqLat;           // Pointer to 1D arrays for trig functions of latitude
    double * cos2Lon, * sin2Lon;            // Pointer to 1D arrays for trig functions of longitude
    double ** latGrad, ** lonGrad;          // Pointer to 2D array solutions
    int i,j;
    int latLen, lonLen;

    factor = pow(omega,2.0)*pow(radius,2.0)*ecc;
    cosM = cos(omega*simulationTime);
    sinM = sin(omega*simulationTime);

    // Get pointer to arrays
    latGrad = DUlat->solution;
    lonGrad = DUlon->solution;

    // Set variables for dUlat solution
    latLen = DUlat->fieldLatLen;
    lonLen = DUlat->fieldLonLen;

    sin2Lat = DUlat->sin2Lat;
    cos2Lon = DUlat->cos2Lon;
    sin2Lon = DUlat->sin2Lon;

    // Solve for dUdlat
    for (i=0; i<latLen; i++) {
        for (j=0; j<lonLen; j++) {
            latGrad[i][j] = -factor*0.75*sin2Lat[i]
                            *(3.*cosM*(1.+cos2Lon[j])
                            + 4*sinM*sin2Lon[j]);
            }
    }

    // Set variable for dUlon solution
    latLen = DUlon->fieldLatLen;
    lonLen = DUlon->fieldLonLen;

    cosSqLat = DUlon->cosSqLat;
    sin2Lon = DUlon->sin2Lon;
    cos2Lon = DUlon->cos2Lon;

    // Solve for dUdlon
    for (i=0; i<latLen; i++) {
        for (j=0; j<lonLen; j++) {
            lonGrad[i][j] = factor* 1.5 * cosSqLat[i]
                            * (4.*sinM * cos2Lon[j]
                            - 3.*cosM * sin2Lon[j]);
        }
    }
};

/* Degree 2 component of the eccentricity-radial tide (see docs)
 *
 * Finds the gradient components of the degree-3 ecc-radial tide for Fields
 * DUlat and DUlon. Only current simulation required is absolute simulation
 * time. The function directly writes the solution arrays in both Field
 * classes. The radial time contains only the P_20 term of eccentricity tide. As
 * the order of potential expansion is zero (m=0), this potential only depends
 * on latitude.
 *
 *      inputs: Field * DUlat           Field for latitude gradient of potential
 *              Field * DUlon           Field for longitude gradient of potential
 *              double simulationTime   Current time in the simulation
 *              double radius           Satellite radius
 *              double omega            Satellite rotational angular speed
 *              double ecc              Satellite eccentricity (must be << 1)
 *
*/
void deg2EccRad(Field * DUlat, Field * DUlon, double simulationTime, double radius, double omega, double ecc) {
    double cosM, factor;                    // cos(Mean anomaly),
    double * sin2Lat;                       // Pointer to 1D arrays for trig functions of latitude
    double ** latGrad, ** lonGrad;          // Pointer to 2D array solutions
    int i,j;
    int latLen, lonLen;
    double val;

    factor = pow(omega,2.0)*pow(radius,2.0)*ecc;
    cosM = cos(omega*simulationTime);

    // Get pointer to arrays
    latGrad = DUlat->solution;
    lonGrad = DUlon->solution;

    // Set variables for dUlat solution
    latLen = DUlat->fieldLatLen;
    lonLen = DUlat->fieldLonLen;

    sin2Lat = DUlat->sin2Lat;

    // Solve for dUdlat
    for (i=0; i<latLen; i++) {
        val = -factor * (2.25 * sin2Lat[i] * cosM);
        for (j=0; j<lonLen; j++) {
            latGrad[i][j] = val;
        }
    }

    // No need to solve for dUdLon as it is always zero in this case

};

/* Degree 2 component of the eccentricity-libration tide (Tyler, 2011)
 *
 * Finds the gradient components of the degree-3 ecc-lib tide for Fields
 * DUlat and DUlon. Only current simulation required is absolute simulation
 * time. The function directly writes the solution arrays in both Field
 * classes. The eccentricity-libration tide contains only the P_22 term from
 * the eccentricity tide.
 *
 *      inputs: Field * DUlat           Field for latitude gradient of potential
 *              Field * DUlon           Field for longitude gradient of potential
 *              double simulationTime   Current time in the simulation
 *              double radius           Satellite radius
 *              double omega            Satellite rotational angular speed
 *              double ecc              Satellite eccentricity (must be << 1)
 *
*/
void deg2EccLib(Field * DUlat, Field * DUlon, double simulationTime, double radius, double omega, double ecc) {
    double cosM, sinM, factor;              // cos(Mean anomaly), sin(Mean anomaly)
    double * sin2Lat, * cosSqLat;           // Pointer to 1D arrays for trig functions of latitude
    double * cos2Lon, * sin2Lon;            // Pointer to 1D arrays for trig functions of longitude
    double ** latGrad, ** lonGrad;          // Pointer to 2D array solutions
    int i,j;
    int latLen, lonLen;

    factor = pow(omega,2.0)*pow(radius,2.0)*ecc;
    cosM = cos(omega*simulationTime);
    sinM = sin(omega*simulationTime);

    // Get pointer to arrays
    latGrad = DUlat->solution;
    lonGrad = DUlon->solution;

    // Set variables for dUlat solution
    latLen = DUlat->fieldLatLen;
    lonLen = DUlat->fieldLonLen;

    sin2Lat = DUlat->sin2Lat;
    cos2Lon = DUlat->cos2Lon;
    sin2Lon = DUlat->sin2Lon;

    // Solve for dUdlat
    for (i=0; i<latLen; i++) {
        for (j=0; j<lonLen; j++) {
            latGrad[i][j] = -factor*0.75*sin2Lat[i]
                            *(3.*cosM*cos2Lon[j]
                            + 4*sinM*sin2Lon[j]);
            }
    }

    // Set variable for dUlon solution
    latLen = DUlon->fieldLatLen;
    lonLen = DUlon->fieldLonLen;

    cosSqLat = DUlon->cosSqLat;
    sin2Lon = DUlon->sin2Lon;
    cos2Lon = DUlon->cos2Lon;

    // Solve for dUdlon
    for (i=0; i<latLen; i++) {
        for (j=0; j<lonLen; j++) {
            lonGrad[i][j] = factor* 1.5 * cosSqLat[i]
                            * (4.*sinM * cos2Lon[j]
                            - 3.*cosM * sin2Lon[j]);
        }
    }
};

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
void deg2Obliq(Field * DUlat, Field * DUlon, double simulationTime, double radius, double omega, double theta)
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

/* Degree 3 component of the obliquity tide (see docs)
 *
 * Finds the gradient components of the degree-3 obliquity tide for Fields
 * DUlat and DUlon. Only current simulation required is absolute simulation
 * time. The function directly writes the solution arrays in both Field
 * classes.
 *
 *      inputs: Field * DUlat           Field for latitude gradient of potential
 *              Field * DUlon           Field for longitude gradient of potential
 *              double simulationTime   Current time in the simulation
 *              double radius           Satellite radius
 *              double smAxis           Satellite semimajor axis
 *              double omega            Satellite rotational angular speed
 *              double theta            Satellite obliquity in radians
 *
*/
void deg3Obliq(Field * DUlat, Field * DUlon, double simulationTime, double radius, double smAxis, double omega, double theta) {
    double cosM, factor;                                // cos(Mean anomaly)
    double * cosLat, * cosSqLat, * sinLat, * sinSqLat;  // Pointer to 1D arrays for trig functions of latitude
    double * cosSqLon, * sin2Lon;                       // Pointer to 1D arrays for trig functions of longitude
    double ** latGrad, ** lonGrad;                      // Pointer to 2D array solutions
    int i,j;
    int latLen, lonLen;

    factor = pow(omega,2.0)*pow(radius,3.0)*theta/smAxis;
    cosM = cos(omega*simulationTime);

    // Get pointer to arrays
    latGrad = DUlat->solution;
    lonGrad = DUlon->solution;

    // Set variables for dUlat solution
    latLen = DUlat->fieldLatLen;
    lonLen = DUlat->fieldLonLen;

    cosLat = DUlat->cosLat;
    sinSqLat = DUlat->sinSqLat;
    cosSqLon = DUlat->cosSqLon;

    // Solve for dUdlat
    for (i=0; i<latLen; i++) {
        for (j=0; j<lonLen; j++) {
            latGrad[i][j] = factor*1.5*cosLat[i]*(5.*cosSqLon[j] * (1.0 - 3.*sinSqLat[i]) - 1.0) * cosM;
        }
    }

    // Set variable for dUlon solution
    latLen = DUlon->fieldLatLen;
    lonLen = DUlon->fieldLonLen;

    cosSqLat = DUlon->cosSqLat;
    sinLat = DUlon->sinLat;
    sin2Lon = DUlon->sin2Lon;

    // Solve for dUdlon
    for (i=0; i<latLen; i++) {
        for (j=0; j<lonLen; j++) {
            lonGrad[i][j] = -factor* 7.5*cosSqLat[i] * sinLat[i] * sin2Lon[j] * cosM;
        }
    }
};
