/* File: tidalPotentials.cpp

   All functions calculate the lateral tidal force (potential gradient) at the
   surface of a spherical body, using analytical expressions. In general:

        inputs: Mesh    grid              mesh object for all spatial info on the grid.
                Array2D<double> soln  Solution array for the calculated forcing
                double simulationTime     Current time in the simulation
                double radius             body radius
                double omega              body rotational angular speed
                double ecc                body eccentricity (must be << 1)
                double obl                body obliquity (in radians)

*/

#include "mesh.h"
#include "globals.h"
#include "tidalPotentials.h"
#include "array2d.h"
#include "spatialOperators.h"
#include <math.h>

// -----------------------------------------------------------------------------
// Full eccentricity tidal forcing, degree-2 (e.g., Matsuyama (2014)) ----------
// -----------------------------------------------------------------------------
void deg2Ecc(Mesh * grid, Array2D<double> & soln, double simulationTime, double radius, double omega, double ecc)
{
    double cosM, sinM, factor;              // cos(Mean anomaly), sin(Mean anomaly)
    double * sin2Lat, * sin2Lon, *cos2Lon, *cosLat;
    double * sinSqLat, * cosSqLat;
    int i,j, node_num;
    double * val;
    double * m;

    node_num = grid->node_num;

    Array1D<double> * scalar_dummy;
    scalar_dummy = new Array1D<double>(node_num);


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
    cos2Lon = &(grid->trig2Lon(0,0));
    sin2Lon = &(grid->trig2Lon(0,1));
    sin2Lat = &(grid->trig2Lat(0,1));

    cosSqLat = &(grid->trigSqLat(0, 0));
    sinSqLat = &(grid->trigSqLat(0, 1));

    // Solve for dUdlon
    for (i=0; i<node_num; i++) {
        j = i*2;

        // soln(i,0) = factor * 1.5 * cosLat[j]   // cosLat here as cos^2(lat)/cos(lat)
        //                 * (4.*sinM * cos2Lon[j]
        //                 - 3.*cosM * sin2Lon[j]);
        //
        // soln(i,1) = -factor*0.75*sin2Lat[j]
        //                 *(3.*cosM*(1.+cos2Lon[j])
        //                 + 4.*sinM*sin2Lon[j]);

        (*scalar_dummy)(i) = 0.75*factor*radius *((1-3*sinSqLat[j])*cosM
                                + cosSqLat[j] * (3*cosM*cos2Lon[j]
                                + 4*sinM*sin2Lon[j]));

    }

    pressureGradient(grid, soln, *scalar_dummy, node_num, -1.0);
    // double a, b;
    // for (i=0; i<node_num; i++) {
    //     j = i*2;
    //
    //     a = factor * 1.5 * cosLat[j]   // cosLat here as cos^2(lat)/cos(lat)
    //                     * (4.*sinM * cos2Lon[j]
    //                     - 3.*cosM * sin2Lon[j]);
    //
    //     b = -factor*0.75*sin2Lat[j]
    //                     *(3.*cosM*(1.+cos2Lon[j])
    //                     + 4.*sinM*sin2Lon[j]);
    //
    //     std::cout<<soln(i,0)<<'\t'<<a<<std::endl;
    //     std::cout<<soln(i,1)<<'\t'<<b<<std::endl;
    //
    //     // (*scalar_dummy)(i) = 0.75*factor*radius *((1-3*sinSqLat[j])*cosM
    //     //                         + cosSqLat[j] * (3*cosM*cos2Lon[j]
    //     //                         + 4*sinM*sin2Lon[j]));
    //
    // }
    delete scalar_dummy;
};


// -----------------------------------------------------------------------------
// Full eccentricity libration tidal forcing, degree-2 (e.g., Matsuyama (2014))
// -----------------------------------------------------------------------------
void deg2EccLib(Mesh * grid, Array2D<double> & soln, double simulationTime, double radius, double omega, double ecc)
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

        soln(i,0) = factor * 1.5 * cosLat[j]   // cosLat here as cos^2(lat)/cos(lat)
                        * (4.*sinM * cos2Lon[j]
                        - 3.*cosM * sin2Lon[j]);

        soln(i,1) = -factor*0.75*sin2Lat[j]
                        *(3.*cosM*cos2Lon[j]
                        + 4.*sinM*sin2Lon[j]);

    }

};


// -----------------------------------------------------------------------------
// Eastward travelling eccentricity tidal forcing, degree-2 (Matsuyama, 2014) --
// -----------------------------------------------------------------------------
void deg2EccEast(Mesh * grid, Array2D<double> & soln, double simulationTime, double radius, double omega, double ecc)
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
        soln(i,0) = factor *5.25 * cosLat[j]   // cosLat here as cos^2(lat)/cos(lat)
                        * (cosM * sin2Lon[j]
                        - sinM * cos2Lon[j]);

        // calculate potential gradient in latitude
        soln(i,1) = factor*2.625*sin2Lat[j]
                        *(cosM*cos2Lon[j]
                        + sinM*sin2Lon[j]);

    }

};


// -----------------------------------------------------------------------------
// Westward travelling eccentricity tidal forcing, degree-2 (Matsuyama, 2014) --
// -----------------------------------------------------------------------------
void deg2EccWest(Mesh * grid, Array2D<double> & soln, double simulationTime, double radius, double omega, double ecc)
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
        soln(i,0) = -factor * 0.75 * cosLat[j]   // cosLat here as cos^2(lat)/cos(lat)
                        * (cosM * sin2Lon[j]
                        + sinM * cos2Lon[j]);

        // calculate potential gradient in latitude
        soln(i,1) = -factor*0.375*sin2Lat[j]
                        *(cosM*cos2Lon[j]
                        - sinM*sin2Lon[j]);

    }

};


// -----------------------------------------------------------------------------
// Radial eccentricity tidal forcing, degree-2 (e.g., Matsuyama 2014, Tyler 2011)
// -----------------------------------------------------------------------------
void deg2EccRad(Mesh * grid, Array2D<double> & soln, double simulationTime, double radius, double omega, double ecc)
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

        soln(i,0) = 0.0;                    // no lon gradient
        soln(i,1) = factor * sin2Lat[j];    // lat gradient
    }

};


// -----------------------------------------------------------------------------
// Full obliquity tidal forcing, degree-2 (e.g., Matsuyama 2014, Tyler 2011)
// -----------------------------------------------------------------------------
void deg2Obliq(Mesh * grid, Array2D<double> & soln, double simulationTime, double radius, double omega, double theta)
{
    double * sinLat, * sinLon, *cosLon, * cos2Lat;
    double * sin2Lat;
    int i,j, node_num;
    double * val;
    double * m;
    double cosM, factor;

    node_num = grid->node_num;

    Array1D<double> * scalar_dummy;
    scalar_dummy = new Array1D<double>(node_num);


    // factor = pow(omega,2.0)*pow(radius,2.0)*ecc;
    factor = -3. * grid->globals->loveReduct.Value() * pow(omega,2.0)*radius*theta;
    cosM = cos(omega*simulationTime);
    // sinM = sin(omega*simulationTime);


    if (grid->globals->surface_type == LID_LOVE || grid->globals->surface_type == LID_MEMBR)
    {
     factor /= radius;
     factor *= pow(radius + grid->globals->shell_thickness.Value(), 2.0)/radius;
    }

    // Assign pointers to start of trig node arrays
    cos2Lat = &(grid->trig2Lat(0,0));
    cosLon = &(grid->trigLon(0,0));
    sinLon = &(grid->trigLon(0,1));
    sinLat = &(grid->trigLat(0,1));
    sin2Lat = &(grid->trig2Lat(0,1));

    // Solve for dUdlon
    for (i=0; i<node_num; i++) {
        j = i*2;
        // calculate potential gradient in longitude
        // soln(i,0) += -factor * sinLat[j] * sinLon[j] * cosM;
        //
        //
        // // calculate potential gradient in latitude
        // soln(i,1) += factor * cos2Lat[j] * cosLon[j] * cosM;

        // Calculate tidal potential
        (*scalar_dummy)(i) = radius*factor/2. * cosM * sin2Lat[j] * cosLon[j];

    }

    pressureGradient(grid, soln, *scalar_dummy, node_num, -1.0);

    delete scalar_dummy;
};


// -----------------------------------------------------------------------------
// Westward travelling obliquity tidal forcing, degree-2 (e.g., Matsuyama 2014)
// -----------------------------------------------------------------------------
void deg2ObliqWest(Mesh * grid, Array2D<double> & soln, double simulationTime, double radius, double omega, double theta)
{
    double * sinLat, * sinLon, *cosLon, * cos2Lat;
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
    cos2Lat = &(grid->trig2Lat(0,0));
    cosLon = &(grid->trigLon(0,0));
    sinLon = &(grid->trigLon(0,1));
    sinLat = &(grid->trigLat(0,1));

    // Solve for dUdlon
    for (i=0; i<node_num; i++) {
        j = i*2;
        // calculate potential gradient in longitude
        soln(i,0) = -factor * sinLat[j] * (sinM * cosLon[j] + sinLon[j] * cosM);


        // calculate potential gradient in latitude
        soln(i,1) = factor * cos2Lat[j] * (cosM * cosLon[j] - sinM * sinLon[j]);

    }
};


// -----------------------------------------------------------------------------
// Eastward travelling obliquity tidal forcing, degree-2 (e.g., Matsuyama 2014)
// -----------------------------------------------------------------------------
void deg2ObliqEast(Mesh * grid, Array2D<double> & soln, double simulationTime, double radius, double omega, double theta)
{
    double * sinLat, * sinLon, *cosLon, * cos2Lat;
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
    cos2Lat = &(grid->trig2Lat(0,0));
    cosLon = &(grid->trigLon(0,0));
    sinLon = &(grid->trigLon(0,1));
    sinLat = &(grid->trigLat(0,1));

    for (i=0; i<node_num; i++) {
        j = i*2;

        // calculate potential gradient in longitude
        soln(i,0) = -factor * sinLat[j] * (cosM * sinLon[j] - sinM * cosLon[j]);


        // calculate potential gradient in latitude
        soln(i,1) = factor * cos2Lat[j] * (sinM * sinLon[j] + cosM * cosLon[j]);

    }
};


// -----------------------------------------------------------------------------
// Full time-varying tidal forcing (ecc + obliq), degree-2 (e.g., Matsuyama 2014)
// -----------------------------------------------------------------------------
void deg2Full(Mesh * grid, Array2D<double> & soln, double simulationTime, double radius, double omega, double theta, double ecc)
{
    double cosM, sinM, factor_lon, factor_lat;              // cos(Mean anomaly), sin(Mean anomaly)
    double * sinLon, * cosLon;
    double * sinLat, * cosLat;
    double * sin2Lat, * cos2Lat;
    double * sin2Lon, * cos2Lon;

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

    // Solve for dUdlon
    for (i=0; i<node_num; i++) {
        j = i*2;
        // calculate potential gradient in longitude
        soln(i,0) = factor_lon * (3. * theta * sinLat[j] * sinLon[j] * cosM
                        + 1.5 * ecc * cosLat[j]
                        * (4. * sinM * cos2Lon[j]
                        - 3. * cosM * sin2Lon[j]));

        // calculate potential gradient in latitude
        soln(i,1) = -factor_lat * (3. * theta * cos2Lat[j] * cosLon[j] * cosM
                        + 0.75 * ecc * sin2Lat[j] * (3. * cosM * (1 + cos2Lon[j])
                        + 4. * sin2Lon[j] * sinM));
    }
};


// -----------------------------------------------------------------------------
// Static + time-varying eccentricity tidal forcing, degree-2 (e.g., Matsuyama 2014)
// -----------------------------------------------------------------------------
void deg2EccTotal(Mesh * grid, Array2D<double> & soln, double simulationTime, double radius, double omega, double theta, double ecc)
{
    double cosM, sinM, factor_lon, factor_lat;              // cos(Mean anomaly), sin(Mean anomaly)
    double * cosLat;
    double * sin2Lat;
    double * sin2Lon, * cos2Lon;
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
    cosLat = &(grid->trigLat(0,0));

    sin2Lat = &(grid->trig2Lat(0,1));
    cos2Lon = &(grid->trig2Lon(0,0));
    sin2Lon = &(grid->trig2Lon(0,1));

    // Solve for dUdlon
    // for (i=0; i<node_num; i++) {
    //     j = i*2;
    //     // factor_lon = factor_lat/cosLat[j];
    //     // if (i < 2) factor_lon = 0.0;
    //     // calculate potential gradient in longitude
    //     soln(i,0) = factor_lon * (-3. * cosLat[j] * sin2Lon[j]
    //                                   + 12. * ecc * cosLat[j] * cos2Lon[j] * sinM);
    //
    //     // calculate potential gradient in latitude
    //     soln(i,1) = factor_lat * (-3. * sin2Lat[j] * cosSqLon[j]
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

        soln(i,0) = factor * 1.5 * cosLat[j]   // cosLat here as cos^2(lat)/cos(lat)
                        * (4.*sinM * cos2Lon[j]
                        - 3.*cosM * sin2Lon[j]);

        soln(i,0) += factor * 1.5 * cosLat[j]* sin2Lon[j];

        soln(i,1) = -factor*0.75*sin2Lat[j]
                        *(3.*cosM*(1.+cos2Lon[j])
                        + 4.*sinM*sin2Lon[j]);

        soln(i,1) += factor * 0.75 * sin2Lat[j] * ( 1. + cos2Lon[j]);

    }

    // avgAtPoles(grid, soln);
};
