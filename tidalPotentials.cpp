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

/* Degree 3 component of the eccentricity tide (see docs)
 *
 * Finds the gradient components of the degree-3 eccentricity tide for Fields
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
 *              double ecc              Satellite eccentricity
 *
*/
void deg3Ecc(Field * DUlat, Field * DUlon, double simulationTime, double radius, double smAxis, double omega, double ecc) {
    double cosM, sinM, factor;                          // cos(Mean anomaly), sin(Mean anomaly)
    double * sinLat, * cosLat, * cosSqLat, * sinSqLat;  // Pointer to 1D arrays for trig functions of latitude
    double * sinLon, * cosLon, * cosSqLon, * cosCubLon; // Pointer to 1D arrays for trig functions of longitude
    double ** latGrad, ** lonGrad;                      // Pointer to 2D array solutions
    int i,j;
    int latLen, lonLen;

    factor = pow(omega,2.0)*pow(radius,3.0)*ecc/smAxis;
    cosM = cos(omega*simulationTime);
    sinM = sin(omega*simulationTime);

    // Get pointer to arrays
    latGrad = DUlat->solution;
    lonGrad = DUlon->solution;

    // Set variables for dUlat solution
    latLen = DUlat->fieldLatLen;
    lonLen = DUlat->fieldLonLen;

    sinLat = DUlat->sinLat;
    sinLon = DUlat->sinLon;
    cosLon = DUlat->cosLon;
    cosSqLat = DUlat->cosSqLat;
    cosSqLon = DUlat->cosSqLon;
    cosCubLon = DUlat->cosCubLon;

    // Solve for dUdlat
    for (i=0; i<latLen; i++) {
        for (j=0; j<lonLen; j++) {
            latGrad[i][j] = -factor * sinLat[i]
                            * (3.*sinM*sinLon[j]
                            * (15.*cosSqLat[i]*cosSqLon[j] - 1.0)
                            + 6. * cosM
                            * (5.*cosSqLat[i]*cosCubLon[j] - cosLon[j]));
        }
    }

    // Set variable for dUlon solution
    latLen = DUlon->fieldLatLen;
    lonLen = DUlon->fieldLonLen;

    cosLat = DUlon->cosLat;
    cosLon = DUlon->cosLon;
    cosSqLat = DUlon->cosSqLat;
    sinSqLat = DUlon->sinSqLat;
    sinLon = DUlon->sinLon;
    cosSqLon = DUlon->cosSqLon;

    // Solve for dUdlon
    for (i=0; i<latLen; i++) {
        for (j=0; j<lonLen; j++) {
            lonGrad[i][j] = -factor * cosLat[i]
                            * (3.*sinM*cosLon[j]
                            * (5.*cosSqLat[i]*(1. - 3.*sinSqLat[i]) - 1.0)
                            + 6. * cosM * sinLon[j]
                            * (1. - 5. * cosSqLat[i] * cosSqLon[j]));
        }
    }
};

// inline void Solver::UpdateFullPotential(void) {
//     double lon = 0;
//     double B = consts->angVel.Value()*simulationTime;
//     double A = 0.25 * pow(consts->angVel.Value(),2)*pow(consts->radius.Value(),2);
//     double sin2Lat = 0;
//     double cos2Lat = 0;
//     double cosLat = 0;
//
//     for (int i = 0; i < dUlat->fieldLatLen; i++) {
//         sin2Lat = dUlat->sin2Lat[i];
//         cos2Lat = dUlat->cos2Lat[i];
//         for (int j = 0; j < dUlat->fieldLonLen; j++) {
//             lon = dUlat->lon[j];
//             dUlatArray[i][j] = consts->theta.Value()* -6 * cos2Lat*(cos(lon - B) + cos(lon + B));
//             dUlatArray[i][j] += consts->e.Value()*(-9.*sin2Lat*cos(B));
//             dUlatArray[i][j] += consts->e.Value()*(-1.5*sin2Lat * (7 * cos(2 * lon - B) - cos(2 * lon + B)));
//             dUlatArray[i][j] *= A;
//         }
//     }
//
//     for (int i = 0; i < dUlon->fieldLatLen; i++) {
//         sin2Lat = dUlon->sin2Lat[i];
//         cosLat = dUlon->cosLat[i];
//         for (int j = 0; j < dUlon->fieldLonLen; j++) {
//             lon = dUlon->lon[j];
//             dUlonArray[i][j] = consts->theta.Value() * 3 * sin2Lat*(sin(lon - B) + sin(lon + B));
//             dUlonArray[i][j] += consts->e.Value() * 3 * cosLat*cosLat*(-7 * sin(2 * lon - B) + sin(2 * lon + B));
//             dUlonArray[i][j] *= A;
//         }
//     }
// }
//
// inline void Solver::UpdateTotalPotential(void) {
//     // Total time-dependent and static parts of the tidal potential, accurate to second order
//     // in both eccentricity and obliquity
//     double factor;
//     // double factor2;
//     double cosLat = 0;
//     double cos2Lat = 0;
//     double sin2Lat = 0;
//     double sinLon = 0;
//     double cosLon = 0;
//     double sin2Lon = 0;
//     double cos2Lon = 0;
//     double lon = 0;
//     double t_omega = consts->angVel.Value()*simulationTime;
//     double e = consts->e.Value();
//     double theta = consts->theta.Value();
//
//     double cos_t_omega = cos(t_omega);
//     double sin_t_omega = sin(t_omega);
//     double cos_2t_omega = cos(2*t_omega);
//     double cos_3t_omega = cos(3*t_omega);
//     double cos_4t_omega = cos(4*t_omega);
//
//     double A = 0;
//     double B = 0;
//     double C = 0;
//     double D = 0;
//     double E = 0;
//     double F = 0;
//     double G = 0;
//
//     factor = 0.03125 * (consts->angVel.Value() * consts->angVel.Value() * consts->radius.Value() * consts->radius.Value());
//     // factor2 = -0.5 * (consts->angVel.Value() * consts->angVel.Value() * consts->radius.Value() * consts->radius.Value());
//
//     A = -(4 + 15 * e*e + 20. * e * cos_t_omega + 43 * e*e * cos_2t_omega);
//     B = 2 * e * (4 + 25 * e * cos_t_omega) * sin_t_omega;
//     C = -2 * theta*theta * (2 + 3*e*e + 6*e*cos_t_omega + 9 * e * e *cos_2t_omega);
//     D = -2 * (2 - 5*e*e + 6*e*cos_t_omega + 17 * e*e*cos_2t_omega);
//     E = 4 * e * (4 + 17*e*cos_t_omega)*sin_t_omega;
//     F = -(-2+theta*theta);
//
//
//     for (int i = 0; i < dUlon->fieldLatLen; i++) {
//         cosLat = dUlon->cosLat[i];
//         sin2Lat = dUlon->sin2Lat[i];
//         for (int j = 0; j < dUlon->fieldLonLen; j++) {
//             sinLon = dUlon->sinLon[j];
//             cosLon = dUlon->cosLon[j];
//             sin2Lon = dUlon->sin2Lon[j];
//             cos2Lon = dUlon->cos2Lon[j];
//             lon = dUlon->lon[j];
//             dUlonArray[i][j] = 12 * theta * sin2Lat * sin_t_omega * (A * sinLon + B * cosLon);
//             dUlonArray[i][j] += 6*cosLat*cosLat * (C * sin(2 * (t_omega + lon)) + F * (D * sin2Lon + E *cos2Lon));
//             dUlonArray[i][j] *= factor;
//
//             //Static part
//             // dUlonArray[i][j] -= factor2 * (3*cosLat*cosLat*sin2Lon);
//         }
//     }
//
//     G = -(2 + 3 * e*e)*(-2 + 3*theta*theta) + 3*e*(4-7*theta*theta)*cos_t_omega;
//     G += 6*(theta*theta + e*e*(3 - 7*theta*theta))*cos_2t_omega;
//     G += 3*e*theta*theta * (7*cos_3t_omega + 17*e*cos_4t_omega);
//     G *= -6.0;
//
//     C *= -0.5;
//     D *= -0.5;
//     E *= 0.5;
//     // F *= 0.5;
//
//
//     for (int i = 0; i < dUlat->fieldLatLen; i++) {
//         cosLat = dUlat->cosLat[i];
//         sin2Lat = dUlat->sin2Lat[i];
//         for (int j = 0; j < dUlat->fieldLonLen; j++) {
//             sinLon = dUlat->sinLon[j];
//             cosLon = dUlat->cosLon[j];
//             sin2Lon = dUlat->sin2Lon[j];
//             cos2Lon = dUlat->cos2Lon[j];
//             lon = dUlat->lon[j];
//             dUlatArray[i][j] = G * sin2Lat;
//             dUlatArray[i][j] += 24 * theta * cos2Lat * sin_t_omega * (-A * cosLon + B * sinLon);
//             dUlatArray[i][j] += -6*sin2Lat * (C * cos(2 * (t_omega + lon)) + F * (D * cos2Lon + E *sin2Lon));
//             dUlatArray[i][j] *= factor;
//
//             //Static part
//             // dUlatArray[i][j] -= factor2 * (3*sin2Lat*cosLon*cosLon);
//         }
//     }
// }
