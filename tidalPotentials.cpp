#include "mesh.h"
#include "globals.h"
#include "tidalPotentials.h"
#include "array2d.h"
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
void deg2Ecc(Mesh * grid, Array2D<double> & velocity, double simulationTime, double radius, double omega, double ecc)
{
    double cosM, sinM, factor;              // cos(Mean anomaly), sin(Mean anomaly)
    double * sin2Lat, * sin2Lon, *cos2Lon, *cosLat;
    int i,j, node_num;
    double * val;
    double * m;
    double k2, h2;

    node_num = grid->node_num;

    k2 = grid->globals->loveK2.Value();
    h2 = grid->globals->loveH2.Value();
    // factor = pow(omega,2.0)*pow(radius,2.0)*ecc;

    factor = (1.0 + k2 - h2) * pow(omega,2.0)*radius*ecc;
    cosM = cos(omega*simulationTime);
    sinM = sin(omega*simulationTime);


    // Assign pointers to start of trig node arrays
    cosLat = &(grid->trigLat(0,0));
    cos2Lon = &(grid->trig2Lon(0,0));
    sin2Lon = &(grid->trig2Lon(0,1));
    sin2Lat = &(grid->trig2Lat(0,1));

    // Solve for dUdlon
    for (i=0; i<node_num; i++) {
        j = i*2;

        velocity(i,0) = factor * 1.5 * cosLat[j]   // cosLat here as cos^2(lat)/cos(lat)
                        * (4.*sinM * cos2Lon[j]
                        - 3.*cosM * sin2Lon[j]);

        velocity(i,1) = -factor*0.75*sin2Lat[j]
                        *(3.*cosM*(1.+cos2Lon[j])
                        + 4.*sinM*sin2Lon[j]);

    }

};

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
void deg2EccEast(Mesh * grid, Array2D<double> & velocity, double simulationTime, double radius, double omega, double ecc)
{
    double cosM, sinM, factor;              // cos(Mean anomaly), sin(Mean anomaly)
    double * sin2Lat, * sin2Lon, *cos2Lon, *cosLat;
    int i,j, node_num;
    double * val;
    double * m;
    double k2, h2;

    node_num = grid->node_num;

    k2 = grid->globals->loveK2.Value();
    h2 = grid->globals->loveH2.Value();
    // factor = pow(omega,2.0)*pow(radius,2.0)*ecc;

    factor = (1.0 + k2 - h2) * pow(omega,2.0)*radius*ecc;
    cosM = cos(omega*simulationTime);
    sinM = sin(omega*simulationTime);


    // Assign pointers to start of trig node arrays
    cosLat = &(grid->trigLat(0,0));
    sin2Lon = &(grid->trig2Lon(0,1));
    cos2Lon = &(grid->trig2Lon(0,0));
    sin2Lat = &(grid->trig2Lat(0,1));

    // Solve for dUdlon
    for (i=0; i<node_num; i++) {
        j = i*2;
        // calculate potential gradient in longitude
        velocity(i,0) = factor *5.25 * cosLat[j]   // cosLat here as cos^2(lat)/cos(lat)
                        * (cosM * sin2Lon[j]
                        - sinM * cos2Lon[j]);

        // calculate potential gradient in latitude
        velocity(i,1) = factor*2.625*sin2Lat[j]
                        *(cosM*cos2Lon[j]
                        + sinM*sin2Lon[j]);

    }

};

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
void deg2EccWest(Mesh * grid, Array2D<double> & velocity, double simulationTime, double radius, double omega, double ecc)
{
    double cosM, sinM, factor;              // cos(Mean anomaly), sin(Mean anomaly)
    double * sin2Lat, * sin2Lon, *cos2Lon, *cosLat;
    int i,j, node_num;
    double * val;
    double * m;
    double k2, h2;

    node_num = grid->node_num;

    k2 = grid->globals->loveK2.Value();
    h2 = grid->globals->loveH2.Value();
    // factor = pow(omega,2.0)*pow(radius,2.0)*ecc;

    factor = (1.0 + k2 - h2) * pow(omega,2.0)*radius*ecc;
    cosM = cos(omega*simulationTime);
    sinM = sin(omega*simulationTime);


    // Assign pointers to start of trig node arrays
    cosLat = &(grid->trigLat(0,0));
    cos2Lon = &(grid->trig2Lon(0,0));
    sin2Lon = &(grid->trig2Lon(0,1));
    sin2Lat = &(grid->trig2Lat(0,1));

    // Solve for dUdlon
    for (i=0; i<node_num; i++) {
        j = i*2;
        // calculate potential gradient in longitude
        velocity(i,0) = -factor * 0.75 * cosLat[j]   // cosLat here as cos^2(lat)/cos(lat)
                        * (cosM * sin2Lon[j]
                        + sinM * cos2Lon[j]);

        // calculate potential gradient in latitude
        velocity(i,1) = -factor*0.375*sin2Lat[j]
                        *(cosM*cos2Lon[j]
                        - sinM*sin2Lon[j]);

    }

};

// /* Degree 2 component of the eccentricity-radial tide (see docs)
//  *
//  * Finds the gradient components of the degree-3 ecc-radial tide for Fields
//  * DUlat and DUlon. Only current simulation required is absolute simulation
//  * time. The function directly writes the solution arrays in both Field
//  * classes. The radial time contains only the P_20 term of eccentricity tide. As
//  * the order of potential expansion is zero (m=0), this potential only depends
//  * on latitude.
//  *
//  *      inputs: Field * DUlat           Field for latitude gradient of potential
//  *              Field * DUlon           Field for longitude gradient of potential
//  *              double simulationTime   Current time in the simulation
//  *              double radius           Satellite radius
//  *              double omega            Satellite rotational angular speed
//  *              double ecc              Satellite eccentricity (must be << 1)
//  *
// */
void deg2EccRad(Mesh * grid, Array2D<double> & velocity, double simulationTime, double radius, double omega, double ecc)
{
    double cosM, factor;              // cos(Mean anomaly), sin(Mean anomaly)
    double * sin2Lat;
    int i, j, node_num;
    double k2, h2;

    k2 = grid->globals->loveK2.Value();
    h2 = grid->globals->loveH2.Value();

    node_num = grid->node_num;

    cosM = cos(omega*simulationTime);
    factor = -2.25 * (1.0 + k2 - h2) *  cosM * pow(omega, 2.0) * radius * ecc;

    sin2Lat = &(grid->trig2Lat(0,1));

    for (i=0; i<node_num; i++) {
        j = i*2;

        velocity(i,0) = 0.0;                    // no lon gradient
        velocity(i,1) = factor * sin2Lat[j];    // lat gradient
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
void deg2Obliq(Mesh * grid, Array2D<double> & velocity, double simulationTime, double radius, double omega, double theta)
{
    double * sinLat, * sinLon, *cosLon, *cosLat, * cos2Lat;
    int i,j, node_num;
    double * val;
    double * m;
    double cosM, factor;
    double k2, h2;

    k2 = grid->globals->loveK2.Value();
    h2 = grid->globals->loveH2.Value();

    node_num = grid->node_num;

    // factor = pow(omega,2.0)*pow(radius,2.0)*ecc;
    factor = 3. * (1.0 + k2 - h2) * pow(omega,2.0)*radius*theta;
    cosM = cos(omega*simulationTime);
    // sinM = sin(omega*simulationTime);

    // Assign pointers to start of trig node arrays
    cosLat = &(grid->trigLat(0,0));
    cos2Lat = &(grid->trig2Lat(0,0));
    cosLon = &(grid->trigLon(0,0));
    sinLon = &(grid->trigLon(0,1));
    sinLat = &(grid->trigLat(0,1));

    // Solve for dUdlon
    for (i=0; i<node_num; i++) {
        j = i*2;
        // calculate potential gradient in longitude
        velocity(i,0) = -factor * sinLat[j] * sinLon[j] * cosM;


        // calculate potential gradient in latitude
        velocity(i,1) = factor * cos2Lat[j] * cosLon[j] * cosM;

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
void deg2ObliqWest(Mesh * grid, Array2D<double> & velocity, double simulationTime, double radius, double omega, double theta)
{
    double * sinLat, * sinLon, *cosLon, *cosLat, * cos2Lat;
    int i,j, node_num;
    double * val;
    double * m;
    double cosM, sinM, factor;
    double k2, h2;

    k2 = grid->globals->loveK2.Value();
    h2 = grid->globals->loveH2.Value();

    node_num = grid->node_num;

    // factor = pow(omega,2.0)*pow(radius,2.0)*ecc;
    factor = 1.5 * (1.0 + k2 - h2) * pow(omega,2.0)*radius*theta;
    cosM = cos(omega*simulationTime);
    sinM = sin(omega*simulationTime);

    // Assign pointers to start of trig node arrays
    cosLat = &(grid->trigLat(0,0));
    cos2Lat = &(grid->trig2Lat(0,0));
    cosLon = &(grid->trigLon(0,0));
    sinLon = &(grid->trigLon(0,1));
    sinLat = &(grid->trigLat(0,1));

    // Solve for dUdlon
    for (i=0; i<node_num; i++) {
        j = i*2;
        // calculate potential gradient in longitude
        velocity(i,0) = -factor * sinLat[j] * (sinM * cosLon[j] + sinLon[j] * cosM);


        // calculate potential gradient in latitude
        velocity(i,1) = factor * cos2Lat[j] * (cosM * cosLon[j] - sinM * sinLon[j]);

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
void deg2ObliqEast(Mesh * grid, Array2D<double> & velocity, double simulationTime, double radius, double omega, double theta)
{
    double * sinLat, * sinLon, *cosLon, *cosLat, * cos2Lat;
    int i,j, node_num;
    double * val;
    double * m;
    double cosM, sinM, factor;
    double k2, h2;

    k2 = grid->globals->loveK2.Value();
    h2 = grid->globals->loveH2.Value();

    node_num = grid->node_num;

    // factor = pow(omega,2.0)*pow(radius,2.0)*ecc;
    factor = 1.5 * (1.0 + k2 - h2) * pow(omega,2.0)*radius*theta;
    cosM = cos(omega*simulationTime);
    sinM = sin(omega*simulationTime);

    // Assign pointers to start of trig node arrays
    cosLat = &(grid->trigLat(0,0));
    cos2Lat = &(grid->trig2Lat(0,0));
    cosLon = &(grid->trigLon(0,0));
    sinLon = &(grid->trigLon(0,1));
    sinLat = &(grid->trigLat(0,1));

    for (i=0; i<node_num; i++) {
        j = i*2;

        // calculate potential gradient in longitude
        velocity(i,0) = -factor * sinLat[j] * (cosM * sinLon[j] - sinM * cosLon[j]);


        // calculate potential gradient in latitude
        velocity(i,1) = factor * cos2Lat[j] * (sinM * sinLon[j] + cosM * cosLon[j]);

    }
};

void deg2Full(Mesh * grid, Array2D<double> & velocity, double simulationTime, double radius, double omega, double theta, double ecc)
{
    double cosM, sinM, factor_lon, factor_lat;              // cos(Mean anomaly), sin(Mean anomaly)
    double * sinLon, * cosLon;
    double * sinLat, * cosLat;
    double * sin2Lat, * cos2Lat;
    double * sin2Lon, * cos2Lon;
    double * cosSqLat;
    double * sinSqLon, * cosSqLon;

    int i,j, node_num;
    double * val;
    double * m;
    double k2, h2;

    node_num = grid->node_num;

    k2 = grid->globals->loveK2.Value();
    h2 = grid->globals->loveH2.Value();
    // factor = pow(omega,2.0)*pow(radius,2.0)*ecc;

    factor_lon = (1.0 + k2 - h2) * pow(omega,2.0)*radius;
    factor_lat = factor_lon;

    cosM = cos(omega*simulationTime);
    sinM = sin(omega*simulationTime);

    // Assign pointers to start of trig node arrays
    cosLon = &(grid->trigLon(0,0));
    sinLon = &(grid->trigLon(0,1));
    cosLat = &(grid->trigLat(0,0));
    sinLat = &(grid->trigLat(0,1));

    cos2Lat = &(grid->trig2Lat(0,0));
    sin2Lat = &(grid->trig2Lat(0,1));
    cos2Lon = &(grid->trig2Lon(0,0));
    sin2Lon = &(grid->trig2Lon(0,1));

    cosSqLat = &(grid->trigSqLat(0,0));
    cosSqLon = &(grid->trigSqLon(0,0));
    sinSqLon = &(grid->trigSqLon(0,1));

    // Solve for dUdlon
    for (i=0; i<node_num; i++) {
        j = i*2;
        // calculate potential gradient in longitude
        velocity(i,0) = factor_lon * (3. * theta * sinLat[j] * sinLon[j] * cosM
                        + 1.5 * ecc * cosLat[j]
                        * (4. * sinM * cos2Lon[j]
                        - 3. * cosM * sin2Lon[j]));

        // calculate potential gradient in latitude
        velocity(i,1) = -factor_lat * (3. * theta * cos2Lat[j] * cosLon[j] * cosM
                        + 0.75 * ecc * sin2Lat[j] * (3. * cosM * (1 + cos2Lon[j])
                        + 4. * sin2Lon[j] * sinM));
    }
};

// OLD VERSION ----------------------------------------------------------------


// /* Degree 2 component of the eccentricity-radial tide (see docs)
//  *
//  * Finds the gradient components of the degree-3 ecc-radial tide for Fields
//  * DUlat and DUlon. Only current simulation required is absolute simulation
//  * time. The function directly writes the solution arrays in both Field
//  * classes. The radial time contains only the P_20 term of eccentricity tide. As
//  * the order of potential expansion is zero (m=0), this potential only depends
//  * on latitude.
//  *
//  *      inputs: Field * DUlat           Field for latitude gradient of potential
//  *              Field * DUlon           Field for longitude gradient of potential
//  *              double simulationTime   Current time in the simulation
//  *              double radius           Satellite radius
//  *              double omega            Satellite rotational angular speed
//  *              double ecc              Satellite eccentricity (must be << 1)
//  *
// */
// void deg2EccRad(Field * DUlat, Field * DUlon, double simulationTime, double radius, double omega, double ecc) {
//     double cosM, factor;                    // cos(Mean anomaly),
//     double * sin2Lat;                       // Pointer to 1D arrays for trig functions of latitude
//     double ** latGrad, ** lonGrad;          // Pointer to 2D array solutions
//     int i,j;
//     int latLen, lonLen;
//     double val;
//
//     factor = pow(omega,2.0)*pow(radius,2.0)*ecc;
//     cosM = cos(omega*simulationTime);
//
//     // Get pointer to arrays
//     latGrad = DUlat->solution;
//     lonGrad = DUlon->solution;
//
//     // Set variables for dUlat solution
//     latLen = DUlat->fieldLatLen;
//     lonLen = DUlat->fieldLonLen;
//
//     sin2Lat = DUlat->sin2Lat;
//
//     // Solve for dUdlat
//     for (i=0; i<latLen; i++) {
//         val = -factor * (2.25 * sin2Lat[i] * cosM);
//         for (j=0; j<lonLen; j++) {
//             latGrad[i][j] = val;
//         }
//     }
//
//     // No need to solve for dUdLon as it is always zero in this case
//
// };
//
// /* Degree 2 component of the eccentricity-libration tide (Tyler, 2011)
//  *
//  * Finds the gradient components of the degree-3 ecc-lib tide for Fields
//  * DUlat and DUlon. Only current simulation required is absolute simulation
//  * time. The function directly writes the solution arrays in both Field
//  * classes. The eccentricity-libration tide contains only the P_22 term from
//  * the eccentricity tide.
//  *
//  *      inputs: Field * DUlat           Field for latitude gradient of potential
//  *              Field * DUlon           Field for longitude gradient of potential
//  *              double simulationTime   Current time in the simulation
//  *              double radius           Satellite radius
//  *              double omega            Satellite rotational angular speed
//  *              double ecc              Satellite eccentricity (must be << 1)
//  *
// */
// void deg2EccLib(Field * DUlat, Field * DUlon, double simulationTime, double radius, double omega, double ecc) {
//     double cosM, sinM, factor;              // cos(Mean anomaly), sin(Mean anomaly)
//     double * sin2Lat, * cosSqLat;           // Pointer to 1D arrays for trig functions of latitude
//     double * cos2Lon, * sin2Lon;            // Pointer to 1D arrays for trig functions of longitude
//     double ** latGrad, ** lonGrad;          // Pointer to 2D array solutions
//     int i,j;
//     int latLen, lonLen;
//
//     factor = pow(omega,2.0)*pow(radius,2.0)*ecc;
//     cosM = cos(omega*simulationTime);
//     sinM = sin(omega*simulationTime);
//
//     // Get pointer to arrays
//     latGrad = DUlat->solution;
//     lonGrad = DUlon->solution;
//
//     // Set variables for dUlat solution
//     latLen = DUlat->fieldLatLen;
//     lonLen = DUlat->fieldLonLen;
//
//     sin2Lat = DUlat->sin2Lat;
//     cos2Lon = DUlat->cos2Lon;
//     sin2Lon = DUlat->sin2Lon;
//
//     // Solve for dUdlat
//     for (i=0; i<latLen; i++) {
//         for (j=0; j<lonLen; j++) {
//             latGrad[i][j] = -factor*0.75*sin2Lat[i]
//                             *(3.*cosM*cos2Lon[j]
//                             + 4*sinM*sin2Lon[j]);
//             }
//     }
//
//     // Set variable for dUlon solution
//     latLen = DUlon->fieldLatLen;
//     lonLen = DUlon->fieldLonLen;
//
//     cosSqLat = DUlon->cosSqLat;
//     sin2Lon = DUlon->sin2Lon;
//     cos2Lon = DUlon->cos2Lon;
//
//     // Solve for dUdlon
//     for (i=0; i<latLen; i++) {
//         for (j=0; j<lonLen; j++) {
//             lonGrad[i][j] = factor* 1.5 * cosSqLat[i]
//                             * (4.*sinM * cos2Lon[j]
//                             - 3.*cosM * sin2Lon[j]);
//         }
//     }
// };
//
// /* Degree 2 component of the obliquity tide (see Tyler 2011, Matsuyama 2014)
//  *
//  * Finds the gradient components of the degree-2 obliquity tide for Fields
//  * DUlat and DUlon. Only current simulation required is absolute simulation
//  * time. The function directly writes the solution arrays in both Field
//  * classes.
//  *
//  *      inputs: Field * DUlat           Field for latitude gradient of potential
//  *              Field * DUlon           Field for longitude gradient of potential
//  *              double simulationTime   Current time in the simulation
//  *              double radius           Satellite radius
//  *              double omega            Satellite rotational angular speed
//  *              double theta            Satellite obliquity in radians
//  *
// */
// void deg2Obliq(Field * DUlat, Field * DUlon, double simulationTime, double radius, double omega, double theta)
// {
//     double cosM, factor;                    // cos(Mean anomaly)
//     double * cosLat, * cos2Lat, * sinLat;   // Pointer to 1D arrays for trig functions of latitude
//     double * cosLon, * sinLon;              // Pointer to 1D arrays for trug functions of longitude
//     double ** latGrad, ** lonGrad;          // Pointer to 2D array solutions
//     int i,j;
//     int latLen, lonLen;
//
//     factor = -3.*pow(omega,2.0)*pow(radius,2.0)*theta;
//     cosM = cos(omega*simulationTime);
//
//     // Get pointer to arrays
//     latGrad = DUlat->solution;
//     lonGrad = DUlon->solution;
//
//     // Set variables for dUlat solution
//     latLen = DUlat->fieldLatLen;
//     lonLen = DUlat->fieldLonLen;
//
//     cos2Lat = DUlat->cos2Lat;
//     cosLon = DUlat->cosLon;
//
//     // Solve for dUdlat
//     for (i=0; i<latLen; i++) {
//         for (j=0; j<lonLen; j++) {
//             latGrad[i][j] = factor*cos2Lat[i]*cosLon[j]*cosM;
//         }
//     }
//
//     // Set variable for dUlon solution
//     latLen = DUlon->fieldLatLen;
//     lonLen = DUlon->fieldLonLen;
//
//     cosLat = DUlon->cosLat;
//     sinLat = DUlon->sinLat;
//     sinLon = DUlon->sinLon;
//
//     // Solve for dUdlon
//     for (i=0; i<latLen; i++) {
//         for (j=0; j<lonLen; j++) {
//             lonGrad[i][j] = -factor*cosLat[i]*sinLat[i]*sinLon[j]*cosM;
//         }
//     }
//
// };
//
// /* Degree 3 component of the obliquity tide (see docs)
//  *
//  * Finds the gradient components of the degree-3 obliquity tide for Fields
//  * DUlat and DUlon. Only current simulation required is absolute simulation
//  * time. The function directly writes the solution arrays in both Field
//  * classes.
//  *
//  *      inputs: Field * DUlat           Field for latitude gradient of potential
//  *              Field * DUlon           Field for longitude gradient of potential
//  *              double simulationTime   Current time in the simulation
//  *              double radius           Satellite radius
//  *              double smAxis           Satellite semimajor axis
//  *              double omega            Satellite rotational angular speed
//  *              double theta            Satellite obliquity in radians
//  *
// */
// void deg3Obliq(Field * DUlat, Field * DUlon, double simulationTime, double radius, double smAxis, double omega, double theta) {
//     double cosM, factor;                                // cos(Mean anomaly)
//     double * cosLat, * cosSqLat, * sinLat, * sinSqLat;  // Pointer to 1D arrays for trig functions of latitude
//     double * cosSqLon, * sin2Lon;                       // Pointer to 1D arrays for trig functions of longitude
//     double ** latGrad, ** lonGrad;                      // Pointer to 2D array solutions
//     int i,j;
//     int latLen, lonLen;
//
//     factor = pow(omega,2.0)*pow(radius,3.0)*theta/smAxis;
//     cosM = cos(omega*simulationTime);
//
//     // Get pointer to arrays
//     latGrad = DUlat->solution;
//     lonGrad = DUlon->solution;
//
//     // Set variables for dUlat solution
//     latLen = DUlat->fieldLatLen;
//     lonLen = DUlat->fieldLonLen;
//
//     cosLat = DUlat->cosLat;
//     sinSqLat = DUlat->sinSqLat;
//     cosSqLon = DUlat->cosSqLon;
//
//     // Solve for dUdlat
//     for (i=0; i<latLen; i++) {
//         for (j=0; j<lonLen; j++) {
//             latGrad[i][j] = factor*1.5*cosLat[i]*(5.*cosSqLon[j] * (1.0 - 3.*sinSqLat[i]) - 1.0) * cosM;
//         }
//     }
//
//     // Set variable for dUlon solution
//     latLen = DUlon->fieldLatLen;
//     lonLen = DUlon->fieldLonLen;
//
//     cosSqLat = DUlon->cosSqLat;
//     sinLat = DUlon->sinLat;
//     sin2Lon = DUlon->sin2Lon;
//
//     // Solve for dUdlon
//     for (i=0; i<latLen; i++) {
//         for (j=0; j<lonLen; j++) {
//             lonGrad[i][j] = -factor* 7.5*cosSqLat[i] * sinLat[i] * sin2Lon[j] * cosM;
//         }
//     }
// };
//
// /* Degree 3 component of the eccentricity tide (see docs)
//  *
//  * Finds the gradient components of the degree-3 eccentricity tide for Fields
//  * DUlat and DUlon. Only current simulation required is absolute simulation
//  * time. The function directly writes the solution arrays in both Field
//  * classes.
//  *
//  *      inputs: Field * DUlat           Field for latitude gradient of potential
//  *              Field * DUlon           Field for longitude gradient of potential
//  *              double simulationTime   Current time in the simulation
//  *              double radius           Satellite radius
//  *              double smAxis           Satellite semimajor axis
//  *              double omega            Satellite rotational angular speed
//  *              double ecc              Satellite eccentricity
//  *
// */
// void deg3Ecc(Field * DUlat, Field * DUlon, double simulationTime, double radius, double smAxis, double omega, double ecc) {
//     double cosM, sinM, factor;                          // cos(Mean anomaly), sin(Mean anomaly)
//     double * sinLat, * cosLat, * cosSqLat, * sinSqLat;  // Pointer to 1D arrays for trig functions of latitude
//     double * sinLon, * cosLon, * cosSqLon, * cosCubLon; // Pointer to 1D arrays for trig functions of longitude
//     double ** latGrad, ** lonGrad;                      // Pointer to 2D array solutions
//     int i,j;
//     int latLen, lonLen;
//
//     factor = pow(omega,2.0)*pow(radius,3.0)*ecc/smAxis;
//     cosM = cos(omega*simulationTime);
//     sinM = sin(omega*simulationTime);
//
//     // Get pointer to arrays
//     latGrad = DUlat->solution;
//     lonGrad = DUlon->solution;
//
//     // Set variables for dUlat solution
//     latLen = DUlat->fieldLatLen;
//     lonLen = DUlat->fieldLonLen;
//
//     sinLat = DUlat->sinLat;
//     sinLon = DUlat->sinLon;
//     cosLon = DUlat->cosLon;
//     cosSqLat = DUlat->cosSqLat;
//     cosSqLon = DUlat->cosSqLon;
//     cosCubLon = DUlat->cosCubLon;
//
//     // Solve for dUdlat
//     for (i=0; i<latLen; i++) {
//         for (j=0; j<lonLen; j++) {
//             latGrad[i][j] = -factor * sinLat[i]
//                             * (3.*sinM*sinLon[j]
//                             * (15.*cosSqLat[i]*cosSqLon[j] - 1.0)
//                             + 6. * cosM
//                             * (5.*cosSqLat[i]*cosCubLon[j] - cosLon[j]));
//         }
//     }
//
//     // Set variable for dUlon solution
//     latLen = DUlon->fieldLatLen;
//     lonLen = DUlon->fieldLonLen;
//
//     cosLat = DUlon->cosLat;
//     cosLon = DUlon->cosLon;
//     cosSqLat = DUlon->cosSqLat;
//     sinSqLat = DUlon->sinSqLat;
//     sinLon = DUlon->sinLon;
//     cosSqLon = DUlon->cosSqLon;
//
//     // Solve for dUdlon
//     for (i=0; i<latLen; i++) {
//         for (j=0; j<lonLen; j++) {
//             lonGrad[i][j] = -factor * cosLat[i]
//                             * (3.*sinM*cosLon[j]
//                             * (5.*cosSqLat[i]*(1. - 3.*sinSqLat[i]) - 1.0)
//                             + 6. * cosM * sinLon[j]
//                             * (1. - 5. * cosSqLat[i] * cosSqLon[j]));
//         }
//     }
// };
