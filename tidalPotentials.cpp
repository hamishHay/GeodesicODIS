#include "mesh.h"
#include "globals.h"
#include "tidalPotentials.h"
#include "array2d.h"
#include "spatialOperators.h"
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
    double *sinLat, * sin2Lat, * sin2Lon, *cos2Lon, *cosLat;
    int i,j, node_num;
    double * val;
    double * m;

    node_num = grid->node_num;
    // factor = pow(omega,2.0)*pow(radius,2.0)*ecc;
    // std::cout<<grid->globals->loveReduct.Value()<<std::endl;

    factor = grid->globals->loveReduct.Value() * pow(omega,2.0)*radius*ecc;

    if (grid->globals->surface_type == LID_LOVE || grid->globals->surface_type == LID_MEMBR)
    {
     factor /= radius;
     factor *= pow(radius + grid->globals->shell_thickness.Value(), 2.0)/radius;
    }

    cosM = cos(omega*simulationTime);
    sinM = sin(omega*simulationTime);


    // Assign pointers to start of trig node arrays
    cosLat = &(grid->trigLat(0,0));
    sinLat = &(grid->trigLat(0,1));
    cos2Lon = &(grid->trig2Lon(0,0));
    sin2Lon = &(grid->trig2Lon(0,1));
    sin2Lat = &(grid->trig2Lat(0,1));

    // Solve for dUdlon
    // for (i=0; i<node_num; i++) {
    //     j = i*2;
    //
    //     velocity(i,0) = factor * 1.5 * cosLat[j]   // cosLat here as cos^2(lat)/cos(lat)
    //                     * (4.*sinM * cos2Lon[j]
    //                     - 3.*cosM * sin2Lon[j]);
    //
    //     velocity(i,1) = -factor*0.75*sin2Lat[j]
    //                     *(3.*cosM*(1.+cos2Lon[j])
    //                     + 4.*sinM*sin2Lon[j]);
    //
    // }

    Array1D<double> * potential;
    potential = new Array1D<double>(node_num);
    for (i=0; i<node_num; i++) {
        j = i*2;

        (*potential)(i) = factor * radius * (-3./4. * (3.*sinLat[j]*sinLat[j] - 1.)* cosM    // cosLat here as cos^2(lat)/cos(lat)
                        + 3./4. * cosLat[j] * cosLat[j] * (4. * sinM * sin2Lon[j] + 3. * cosM * cos2Lon[j]));

        velocity(i, 0) = 0.0;
        velocity(i, 1) = 0.0;
    }

    pressureGradient(grid, velocity, *potential, node_num, -1.0);

    delete potential;

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
void deg2EccLib(Mesh * grid, Array2D<double> & velocity, double simulationTime, double radius, double omega, double ecc)
{
    double cosM, sinM, factor;              // cos(Mean anomaly), sin(Mean anomaly)
    double * sin2Lat, * sin2Lon, *cos2Lon, *cosLat;
    int i,j, node_num;
    double * val;
    double * m;

    node_num = grid->node_num;
    // factor = pow(omega,2.0)*pow(radius,2.0)*ecc;

    // std::cout<<grid->globals->loveReduct.Value()<<std::endl;
    factor = grid->globals->loveReduct.Value() * pow(omega,2.0)*radius*ecc;
    cosM = cos(omega*simulationTime);
    sinM = sin(omega*simulationTime);

    if (grid->globals->surface_type == LID_LOVE || grid->globals->surface_type == LID_MEMBR)
    {
     factor /= radius;
     factor *= pow(radius + grid->globals->shell_thickness.Value(), 2.0)/radius;
    }


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
                        *(3.*cosM*cos2Lon[j]
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

    node_num = grid->node_num;
    // factor = pow(omega,2.0)*pow(radius,2.0)*ecc;

    factor = grid->globals->loveReduct.Value() * pow(omega,2.0)*radius*ecc;
    cosM = cos(omega*simulationTime);
    sinM = sin(omega*simulationTime);

    if (grid->globals->surface_type == LID_LOVE || grid->globals->surface_type == LID_MEMBR)
    {
     factor /= radius;
     factor *= pow(radius + grid->globals->shell_thickness.Value(), 2.0)/radius;
    }


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

    node_num = grid->node_num;
    // factor = pow(omega,2.0)*pow(radius,2.0)*ecc;

    factor = grid->globals->loveReduct.Value() * pow(omega,2.0)*radius*ecc;
    cosM = cos(omega*simulationTime);
    sinM = sin(omega*simulationTime);

    if (grid->globals->surface_type == LID_LOVE || grid->globals->surface_type == LID_MEMBR)
    {
     factor /= radius;
     factor *= pow(radius + grid->globals->shell_thickness.Value(), 2.0)/radius;
    }


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

    node_num = grid->node_num;

    cosM = cos(omega*simulationTime);
    factor = -2.25 * grid->globals->loveReduct.Value() *  cosM * pow(omega, 2.0) * radius * ecc;

    if (grid->globals->surface_type == LID_LOVE || grid->globals->surface_type == LID_MEMBR)
    {
     factor /= radius;
     factor *= pow(radius + grid->globals->shell_thickness.Value(), 2.0)/radius;
    }

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
    double * sinLat, * sin2Lat, * sinLon, *cosLon, *cosLat, * cos2Lat;
    int i,j, node_num;
    double * val;
    double * m;
    double cosM, factor;

    node_num = grid->node_num;

    factor = -3. * grid->globals->loveReduct.Value() * pow(omega,2.0)*radius*theta;
    cosM = cos(omega*simulationTime);

    if (grid->globals->surface_type == LID_LOVE || grid->globals->surface_type == LID_MEMBR)
    {
     factor /= radius;
     factor *= pow(radius + grid->globals->shell_thickness.Value(), 2.0)/radius;
    }

    // Assign pointers to start of trig node arrays
    cosLat = &(grid->trigLat(0,0));
    cos2Lat = &(grid->trig2Lat(0,0));
    sin2Lat = &(grid->trig2Lat(0,1));
    cosLon = &(grid->trigLon(0,0));
    sinLon = &(grid->trigLon(0,1));
    sinLat = &(grid->trigLat(0,1));

    // Solve for dUdlon
    // for (i=0; i<node_num; i++) {
    //     j = i*2;
    //     // calculate potential gradient in longitude
    //     velocity(i,0) = -factor * sinLat[j] * sinLon[j] * cosM;
    //
    //
    //     // calculate potential gradient in latitude
    //     velocity(i,1) = factor * cos2Lat[j] * cosLon[j] * cosM;
    //
    // }

    Array1D<double> * potential;
    potential = new Array1D<double>(node_num);
    for (i=0; i<node_num; i++) {
        j = i*2;

        (*potential)(i) = factor/2.0 * radius * sin2Lat[j] * cosLon[j] * cosM;

        velocity(i, 0) = 0.0;
        velocity(i, 1) = 0.0;
    }

    pressureGradient(grid, velocity, *potential, node_num, -1.0);

    delete potential;
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

    node_num = grid->node_num;

    // factor = pow(omega,2.0)*pow(radius,2.0)*ecc;
    factor = 1.5 * grid->globals->loveReduct.Value() * pow(omega,2.0)*radius*theta;
    cosM = cos(omega*simulationTime);
    sinM = sin(omega*simulationTime);

    if (grid->globals->surface_type == LID_LOVE || grid->globals->surface_type == LID_MEMBR)
    {
     factor /= radius;
     factor *= pow(radius + grid->globals->shell_thickness.Value(), 2.0)/radius;
    }

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

    node_num = grid->node_num;

    // factor = pow(omega,2.0)*pow(radius,2.0)*ecc;
    factor = 1.5 * grid->globals->loveReduct.Value() * pow(omega,2.0)*radius*theta;

    if (grid->globals->surface_type == LID_LOVE || grid->globals->surface_type == LID_MEMBR)
    {
     factor /= radius;
     factor *= pow(radius + grid->globals->shell_thickness.Value(), 2.0)/radius;
    }

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

    node_num = grid->node_num;
    // factor = pow(omega,2.0)*pow(radius,2.0)*ecc;

    factor_lon = grid->globals->loveReduct.Value() * pow(omega,2.0)*radius;

    if (grid->globals->surface_type == LID_LOVE || grid->globals->surface_type == LID_MEMBR)
    {
     factor_lon /= radius;
     factor_lon *= pow(radius + grid->globals->shell_thickness.Value(), 2.0)/radius;
    }

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

void deg2EccTotal(Mesh * grid, Array2D<double> & velocity, double simulationTime, double radius, double omega, double theta, double ecc)
{
    double cosM, sinM, factor_lon, factor_lat;              // cos(Mean anomaly), sin(Mean anomaly)
    double * sinLon, * cosLon;
    double * sinLat, * cosLat;
    double * sin2Lat, * cos2Lat;
    double * sin2Lon, * cos2Lon;
    double * cosSqLat;
    double * sinSqLon, * cosSqLon;
    double factor;

    int i,j, node_num;
    double * val;
    double * m;

    node_num = grid->node_num;
    //
    // factor_lon = 0.5 * grid->globals->loveReduct.Value() * pow(omega,2.0)*radius;
    //
    // if (grid->globals->surface_type == LID_LOVE || grid->globals->surface_type == LID_MEMBR)
    // {
    //  factor_lon /= radius;
    //  factor_lon *= pow(radius + grid->globals->shell_thickness.Value(), 2.0)/radius;
    // }

    cosM = cos(omega*simulationTime);
    sinM = sin(omega*simulationTime);
    //
    // factor_lon *= (1. + 3.*ecc*cosM);
    // factor_lat = factor_lon;
    //
    // // Assign pointers to start of trig node arrays
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
    // for (i=0; i<node_num; i++) {
    //     j = i*2;
    //     // factor_lon = factor_lat/cosLat[j];
    //     // if (i < 2) factor_lon = 0.0;
    //     // calculate potential gradient in longitude
    //     velocity(i,0) = factor_lon * (-3. * cosLat[j] * sin2Lon[j]
    //                                   + 12. * ecc * cosLat[j] * cos2Lon[j] * sinM);
    //
    //     // calculate potential gradient in latitude
    //     velocity(i,1) = factor_lat * (-3. * sin2Lat[j] * cosSqLon[j]
    //                                   - 12. * ecc * sin2Lat[j] * sinLon[j] * cosLon[j] * sinM);
    // }

    factor = grid->globals->loveReduct.Value() * pow(omega,2.0)*radius*ecc;

    if (grid->globals->surface_type == LID_LOVE || grid->globals->surface_type == LID_MEMBR)
    {
     factor /= radius;
     factor *= pow(radius + grid->globals->shell_thickness.Value(), 2.0)/radius;
    }

    for (i=0; i<node_num; i++) {
        j = i*2;

        velocity(i,0) = factor * 1.5 * cosLat[j]   // cosLat here as cos^2(lat)/cos(lat)
                        * (4.*sinM * cos2Lon[j]
                        - 3.*cosM * sin2Lon[j]);

        velocity(i,0) += factor * 1.5 * cosLat[j]* sin2Lon[j];

        velocity(i,1) = -factor*0.75*sin2Lat[j]
                        *(3.*cosM*(1.+cos2Lon[j])
                        + 4.*sinM*sin2Lon[j]);

        velocity(i,1) += factor * 0.75 * sin2Lat[j] * ( 1. + cos2Lon[j]);

    }

    // avgAtPoles(grid, velocity);
};
