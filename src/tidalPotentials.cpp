/* File: tidalPotentials.cpp
 *
 * All functions calculate the lateral tidal force (potential gradient) at the
 * surface of a spherical body, using analytical expressions. In general:
 *
 *      inputs: Mesh  * grid              mesh object for all spatial info on the grid.
 *              Array2D<double> soln  Solution array for the calculated forcing
 *              double simulationTime     Current time in the simulation
 *              double radius             body radius
 *              double omega              body rotational angular speed
 *              double ecc                body eccentricity (must be << 1)
 *              double obl                body obliquity (in radians)
 *
*/

// TODO - REWRITE ALL FORCINGS IN TERMS OF THE ABSOLUTE POTENTIAL INSTEAD OF THEIR
// GRADIENTS!!!!!

#include "mesh.h"
#include "globals.h"
#include "tidalPotentials.h"
#include "array2d.h"
#include "spatialOperators.h"
#include <math.h>
#include <iomanip>


void forcing(Globals * consts, Mesh * grid, Array1D<double> & potential, int forcing_type, double time, double ecc, double obl)
{
    int i,j, node_num;
    node_num = grid->node_num;

    double factor, factor2;
    double radius, omega;

    double * sinLon, * cosLon;
    double * sin2Lon, *cos2Lon;
    double * cosLat, * sinLat;
    double * sin2Lat, *cos2Lat;
    double * sinSqLat, * cosSqLat;

    double cosM, sinM;
    double cos2M, sin2M;

    radius = consts->radius.Value();
    omega = consts->angVel.Value();



    if (consts->surface_type == LID_LOVE || consts->surface_type == LID_MEMBR)
    {
     radius += consts->shell_thickness.Value();
    }

    cosM = cos(omega*time);
    sinM = sin(omega*time);
    cos2M = cos(2*omega*time);
    sin2M = sin(2*omega*time);

    double cos3M = cos(3*omega*time);
    double cos4M = cos(4*omega*time);

    // Assign pointers to start of trig node arrays
    cosLon = &(grid->trigLon(0,0));
    sinLon = &(grid->trigLon(0,1));
    cosLat = &(grid->trigLat(0,0));
    sinLat = &(grid->trigLat(0,1));


    cos2Lon = &(grid->trig2Lon(0,0));
    sin2Lon = &(grid->trig2Lon(0,1));
    sin2Lat = &(grid->trig2Lat(0,1));

    cos2Lat = &(grid->trig2Lat(0,0));

    cosSqLat = &(grid->trigSqLat(0, 0));
    sinSqLat = &(grid->trigSqLat(0, 1));


    switch (forcing_type)
    {
      case ECC:
      {
          factor = 0.75 * consts->loveReduct.Value() * pow(omega,2.0)*pow(radius,2.0)*ecc;

          #pragma omp parallel for
          for (i=0; i<node_num; ++i) {
              j = i*2;

              potential(i) = factor*((1. - 3.*sinSqLat[i*2])*cosM
                                      + cosSqLat[i*2] * (3.*cosM*cos2Lon[i*2]
                                      + 4.*sinM*sin2Lon[i*2]));


            double U2 = 9./2. * cos2M * (0.75*cosSqLat[i*2] - 0.5);
            U2 += 17* 3/8. * cosSqLat[i*2] * (sin2M*sin2Lon[i*2] + cos2M*cos2Lon[2*i]);
            U2 *= factor * 4./3. * ecc;
            potential(i) += U2;

          }
      }
      break;

      case OBLIQ:
      {
          factor = -3./2. * consts->loveReduct.Value() * pow(omega,2.0)*pow(radius,2.0)*obl;

          #pragma omp parallel for
          for (i=0; i<node_num; ++i) {
              j = i*2;

              potential(i) = factor * cosM * sin2Lat[j] * cosLon[j];
          }

      }
      break;

      case OBLIQ_WEST:
      {
          factor = 0.5 * consts->loveReduct.Value() * pow(omega,2.0)*pow(radius,2.0)*obl;

          #pragma omp parallel for
          for (i=0; i<node_num; ++i) {
              j = i*2;

              potential(i) = -3*factor * sinLat[j]*cosLat[j]*(cosLon[j]*cosM - sinLon[j]*sinM);
          }

      }
      break;

      case FULL2:
      {
          factor = 1/32. * consts->loveReduct.Value() * pow(omega,2.0)*pow(radius,2.0);
          for (i=0; i<node_num; i++) {
              j = i*2;
              double T1, T2, T3;

              T1 = 3.*ecc*(4. - 7.*obl*obl) * cosM + 6*(obl*obl + ecc*ecc *(3 - 7*obl*obl) ) * cos2M;
              T1 += 3*ecc*obl*obl * ( 7*cos3M + 17*ecc*cos4M );
              T1 *= -(1 - 3* cos2Lat[j]);

              T2 = (4 + 15*ecc*ecc + 20*ecc*cosM + 43*ecc*ecc*cos2M)*cosLon[j];
              T2 += 2*ecc * (4 + 25*ecc*cosM) * sinM * sinLon[j];
              T2 *= 24 * obl * cosLat[j] * sinLat[j] * sinM;

              T3 = obl*obl * (2 + 3*ecc*ecc + 6*ecc*cosM + 9*ecc*ecc*cos2M)*(cosM*cosLon[j] + sinM*sinLon[j]);
              T3 += -(obl*obl - 2) * ( (6*ecc*cosM + 17*ecc*ecc*cos2M)*cos2Lon[j] + 2*ecc*(4+ 17*ecc*cosM) * sinM *sin2Lon[j] );
              T3 *= 6*cosSqLat[j];

              potential(i) = factor*(T1 + T2 + T3);

          }
      }
      break;

      case FULL:
      {
          factor = 0.75 * consts->loveReduct.Value() * pow(omega,2.0)*pow(radius,2.0)*ecc;

          factor2 = -3./2. * consts->loveReduct.Value() * pow(omega,2.0)*pow(radius,2.0)*obl;

          #pragma omp parallel for
          for (i=0; i<node_num; ++i) {
              j = i*2;

              potential(i) = factor*((1-3*sinSqLat[j])*cosM
                                      + cosSqLat[j] * (3*cosM*cos2Lon[j]
                                      + 4*sinM*sin2Lon[j]))
                            + factor2 * cosM * sin2Lat[j] * cosLon[j];
          }
      }
      break;

      case PLANET:
      {
          double a1, a2;  //semimajor axes
          double n1, n2;  //mean motions
          double m2;
          double cosnt, sinnt;
          double cosphi, sinphi;
          double dist;
          double cosgam;

          m2 = 8.931938e+22;  // Io
          // m2 = 1.4815e23;   // Ganymede
          // m2 = 4.799e22;  //Europa

          a1 = 421800000.0; // Io
          // a2 = 1.074e9;   // Ganymede
          a2 = grid->globals->a.Value();
          // a2 = 671100000.0; // Europa

          // n2 = 4.11e-5;   // Io
          // n2 = 1.016e-5;     // Ganymede
          // n2 = 2.05e-5;     // Europa
          n2 = grid->globals->angVel.Value();
          n1 = n2*2.0;

          double nij = (n1-n2);

          cosnt = cos(nij*time);
          sinnt = sin(nij*time);
          double p = pow(a1, 2.0) + pow(a2, 2.0) - 2.*a1*a2*cosnt;
          dist = sqrt(pow(a1, 2.0) + pow(a2, 2.0) - 2.*a1*a2*cosnt);

          cosphi = (a1 - a2*cosnt);
          sinphi = a2*sinnt;

          factor = 0.5 * 6.67408e-11 * m2 * pow(radius/p, 2.0)/sqrt(p);

          // std::cout<<nij<<std::endl;

          // if (i=0)
          #pragma omp parallel for
          for (i=0; i<node_num; ++i) {
              j = i*2;

              cosgam = cosLat[j] * (cosLon[j]*cosphi + sinLon[j]*sinphi);
              potential(i) = factor * (3. * pow(cosgam, 2.0) - p);
              // if (i==0) std::cout<<pow(cosgam, 2.0)<<std::endl;
          }
        }
        break;

        case GENERAL:
        {
            double a1, a2;  //semimajor axes
            double n1, n2;  //mean motions
            double m2;
            // double cosnt, sinnt;
            double cosphi, sinphi;
            double dist;
            double cosgam;
            // double U20, U22a, U22b

            double * P20 = &(grid->Pbar_20(0));
            double * P22 = &(grid->Pbar_22(0));

            double a20, a22, b22;
            int q_max = grid->globals->freq.Value();


            a20 = -9.991062592506567e-05*2.0;
            a22 =   0.0012049514054640862*2.0;
            b22 =  -0.0012017622891820471*2.0;

            n2 = grid->globals->angVel.Value();
            n1 = n2*2.0;

            double nij = (n1-n2);

            #pragma omp parallel for
            for (i=0; i<node_num; i++) {
                j = i*2;
                double U20=0.0;
                double U22a=0.0;
                double U22b=0.0;

                double cos_qnt = cos( (double)q_max * nij * time);
                double sin_qnt = sin( (double)q_max * nij * time);

                U20 += (a20 * cos_qnt);
                U22a += (a22 * cos_qnt);
                U22b += (b22 * sin_qnt);


                U20 = (U20)*P20[i];
                U22a = (U22a)*P22[i]*cos2Lon[j];
                U22b = (U22b)*P22[i]*sin2Lon[j];

                potential(i) = (U20 + U22a + U22b);//*grid->globals->loveReduct.Value();

                // cosgam = cosLat[j] * (cosLon[j]*cosphi + sinLon[j]*sinphi);
                // double test = factor * (3. * pow(cosgam, 2.0) - 1.0)*grid->globals->loveReduct.Value();

                // std::cout<<test<<' '<<(*scalar_dummy)(i)<<' '<<test/(*scalar_dummy)(i)<<std::endl;
            }
      }
      break;
    }

    // switch (globals->tide_type)
    // {
    //     case ECC:
    //         deg2Ecc(grid, dvdt, current_time, r, omega, e);
    //         break;
    //     case ECC_LIB:
    //         deg2EccLib(grid, dvdt, current_time, r, omega, e);
    //         break;
    //     case ECC_WEST:
    //         deg2EccWest(grid, dvdt, current_time, r, omega, e);
    //         break;
    //     case ECC_EAST:
    //         deg2EccEast(grid, dvdt, current_time, r, omega, e);
    //         break;
    //     case ECC_RAD:
    //         deg2EccRad(grid, dvdt, current_time, r, omega, e);
    //         break;
    //     case OBLIQ:
    //         deg2Obliq(grid, dvdt, current_time, r, omega, obliq);
    //         break;
    //     case OBLIQ_WEST:
    //         deg2ObliqWest(grid, dvdt, current_time, r, omega, obliq);
    //         break;
    //     case OBLIQ_EAST:
    //         deg2ObliqEast(grid, dvdt, current_time, r, omega, obliq);
    //         break;
    //     case FULL:
    //         deg2Full(grid, dvdt, current_time, r, omega, obliq, e);
    //         // deg2Ecc(grid, dvdt, current_time, r, omega, e);
    //         // deg2Obliq(grid, dvdt, current_time, r, omega, obliq);
    //         break;
    //     case PLANET:
    //         deg2Planet(grid, dvdt, *forcing_potential, current_time, r);
    //         break;
    //     case PLANET_OBL:
    //         deg2PlanetObl(grid, dvdt, current_time, r, obliq);
    //         break;
    //     case GENERAL:
    //         deg2General(grid, dvdt, current_time, r, p_tm1);
    //         break;
    // }
};


// -----------------------------------------------------------------------------
// Full eccentricity tidal forcing, degree-2 (e.g., Matsuyama (2014)) ----------
// -----------------------------------------------------------------------------
// void deg2Ecc(Mesh * grid, Array2D<double> & soln, double simulationTime, double radius, double omega, double ecc)
// {
//     double cosM, sinM, factor;              // cos(Mean anomaly), sin(Mean anomaly)
//     double * sin2Lat, * sin2Lon, *cos2Lon, *cosLat;
//     double * sinSqLat, * cosSqLat;
//     int i,j, node_num;
//     double * val;
//     double * m;
//
//     node_num = grid->node_num;
//
//
//     Array1D<double> * scalar_dummy;
//     scalar_dummy = new Array1D<double>(node_num);
//
//
//
//     if (grid->globals->surface_type == LID_LOVE || grid->globals->surface_type == LID_MEMBR)
//     {
//      radius += grid->globals->shell_thickness.Value();
//      // factor /= radius;
//      // factor *= pow(radius + grid->globals->shell_thickness.Value(), 2.0)/radius;
//     }
//     factor = grid->globals->loveReduct.Value() * pow(omega,2.0)*radius*ecc;
//
//     cosM = cos(omega*simulationTime);
//     sinM = sin(omega*simulationTime);
//
//     // Assign pointers to start of trig node arrays
//     cosLat = &(grid->trigLat(0,0));
//     cos2Lon = &(grid->trig2Lon(0,0));
//     sin2Lon = &(grid->trig2Lon(0,1));
//     sin2Lat = &(grid->trig2Lat(0,1));
//
//     cosSqLat = &(grid->trigSqLat(0, 0));
//     sinSqLat = &(grid->trigSqLat(0, 1));
//
//     // Solve for dUdlon
//     for (i=0; i<node_num; i++) {
//         j = i*2;
//
//         // soln(i,0) = factor * 1.5 * cosLat[j]   // cosLat here as cos^2(lat)/cos(lat)
//         //                 * (4.*sinM * cos2Lon[j]
//         //                 - 3.*cosM * sin2Lon[j]);
//         //
//         // soln(i,1) = -factor*0.75*sin2Lat[j]
//         //                 *(3.*cosM*(1.+cos2Lon[j])
//         //                 + 4.*sinM*sin2Lon[j]);
//
//         (*scalar_dummy)(i) = 0.75*factor*radius *((1-3*sinSqLat[j])*cosM
//                                 + cosSqLat[j] * (3*cosM*cos2Lon[j]
//                                 + 4*sinM*sin2Lon[j]));
//
//     }
//
//     pressureGradient(grid, soln, *scalar_dummy, node_num, -1.0);
//
//     delete scalar_dummy;
// };
//
//
// // -----------------------------------------------------------------------------
// // Full eccentricity libration tidal forcing, degree-2 (e.g., Matsuyama (2014))
// // -----------------------------------------------------------------------------
// void deg2EccLib(Mesh * grid, Array2D<double> & soln, double simulationTime, double radius, double omega, double ecc)
// {
//     double cosM, sinM, factor;              // cos(Mean anomaly), sin(Mean anomaly)
//     double * sin2Lat, * sin2Lon, *cos2Lon, *cosLat;
//     int i,j, node_num;
//     double * val;
//     double * m;
//
//     node_num = grid->node_num;
//     // factor = pow(omega,2.0)*pow(radius,2.0)*ecc;
//
//     // std::cout<<grid->globals->loveReduct.Value()<<std::endl;
//     factor = grid->globals->loveReduct.Value() * pow(omega,2.0)*radius*ecc;
//     cosM = cos(omega*simulationTime);
//     sinM = sin(omega*simulationTime);
//
//     if (grid->globals->surface_type == LID_LOVE || grid->globals->surface_type == LID_MEMBR)
//     {
//      factor /= radius;
//      factor *= pow(radius + grid->globals->shell_thickness.Value(), 2.0)/radius;
//     }
//
//
//     // Assign pointers to start of trig node arrays
//     cosLat = &(grid->trigLat(0,0));
//     cos2Lon = &(grid->trig2Lon(0,0));
//     sin2Lon = &(grid->trig2Lon(0,1));
//     sin2Lat = &(grid->trig2Lat(0,1));
//
//     // Solve for dUdlon
//     for (i=0; i<node_num; i++) {
//         j = i*2;
//
//         soln(i,0) = factor * 1.5 * cosLat[j]   // cosLat here as cos^2(lat)/cos(lat)
//                         * (4.*sinM * cos2Lon[j]
//                         - 3.*cosM * sin2Lon[j]);
//
//         soln(i,1) = -factor*0.75*sin2Lat[j]
//                         *(3.*cosM*cos2Lon[j]
//                         + 4.*sinM*sin2Lon[j]);
//
//     }
//
// };
//
//
// // -----------------------------------------------------------------------------
// // Eastward travelling eccentricity tidal forcing, degree-2 (Matsuyama, 2014) --
// // -----------------------------------------------------------------------------
// void deg2EccEast(Mesh * grid, Array2D<double> & soln, double simulationTime, double radius, double omega, double ecc)
// {
//     double cosM, sinM, factor;              // cos(Mean anomaly), sin(Mean anomaly)
//     double * sin2Lat, * sin2Lon, *cos2Lon, *cosLat;
//     int i,j, node_num;
//     double * val;
//     double * m;
//
//     node_num = grid->node_num;
//     // factor = pow(omega,2.0)*pow(radius,2.0)*ecc;
//
//     factor = grid->globals->loveReduct.Value() * pow(omega,2.0)*radius*ecc;
//     cosM = cos(omega*simulationTime);
//     sinM = sin(omega*simulationTime);
//
//     if (grid->globals->surface_type == LID_LOVE || grid->globals->surface_type == LID_MEMBR)
//     {
//      factor /= radius;
//      factor *= pow(radius + grid->globals->shell_thickness.Value(), 2.0)/radius;
//     }
//
//
//     // Assign pointers to start of trig node arrays
//     cosLat = &(grid->trigLat(0,0));
//     sin2Lon = &(grid->trig2Lon(0,1));
//     cos2Lon = &(grid->trig2Lon(0,0));
//     sin2Lat = &(grid->trig2Lat(0,1));
//
//     // Solve for dUdlon
//     for (i=0; i<node_num; i++) {
//         j = i*2;
//         // calculate potential gradient in longitude
//         soln(i,0) = factor *5.25 * cosLat[j]   // cosLat here as cos^2(lat)/cos(lat)
//                         * (cosM * sin2Lon[j]
//                         - sinM * cos2Lon[j]);
//
//         // calculate potential gradient in latitude
//         soln(i,1) = factor*2.625*sin2Lat[j]
//                         *(cosM*cos2Lon[j]
//                         + sinM*sin2Lon[j]);
//
//     }
//
// };
//
//
// // -----------------------------------------------------------------------------
// // Westward travelling eccentricity tidal forcing, degree-2 (Matsuyama, 2014) --
// // -----------------------------------------------------------------------------
// void deg2EccWest(Mesh * grid, Array2D<double> & soln, double simulationTime, double radius, double omega, double ecc)
// {
//     double cosM, sinM, factor;              // cos(Mean anomaly), sin(Mean anomaly)
//     double * sin2Lat, * sin2Lon, *cos2Lon, *cosLat;
//     int i,j, node_num;
//     double * val;
//     double * m;
//
//     node_num = grid->node_num;
//     // factor = pow(omega,2.0)*pow(radius,2.0)*ecc;
//
//     factor = grid->globals->loveReduct.Value() * pow(omega,2.0)*radius*ecc;
//     cosM = cos(omega*simulationTime);
//     sinM = sin(omega*simulationTime);
//
//     if (grid->globals->surface_type == LID_LOVE || grid->globals->surface_type == LID_MEMBR)
//     {
//      factor /= radius;
//      factor *= pow(radius + grid->globals->shell_thickness.Value(), 2.0)/radius;
//     }
//
//
//     // Assign pointers to start of trig node arrays
//     cosLat = &(grid->trigLat(0,0));
//     cos2Lon = &(grid->trig2Lon(0,0));
//     sin2Lon = &(grid->trig2Lon(0,1));
//     sin2Lat = &(grid->trig2Lat(0,1));
//
//     // Solve for dUdlon
//     for (i=0; i<node_num; i++) {
//         j = i*2;
//         // calculate potential gradient in longitude
//         soln(i,0) = -factor * 0.75 * cosLat[j]   // cosLat here as cos^2(lat)/cos(lat)
//                         * (cosM * sin2Lon[j]
//                         + sinM * cos2Lon[j]);
//
//         // calculate potential gradient in latitude
//         soln(i,1) = -factor*0.375*sin2Lat[j]
//                         *(cosM*cos2Lon[j]
//                         - sinM*sin2Lon[j]);
//
//     }
//
// };
//
//
// // -----------------------------------------------------------------------------
// // Radial eccentricity tidal forcing, degree-2 (e.g., Matsuyama 2014, Tyler 2011)
// // -----------------------------------------------------------------------------
// void deg2EccRad(Mesh * grid, Array2D<double> & soln, double simulationTime, double radius, double omega, double ecc)
// {
//     double cosM, factor;              // cos(Mean anomaly), sin(Mean anomaly)
//     double * sin2Lat;
//     int i, j, node_num;
//
//     node_num = grid->node_num;
//
//     cosM = cos(omega*simulationTime);
//     factor = -2.25 * grid->globals->loveReduct.Value() *  cosM * pow(omega, 2.0) * radius * ecc;
//
//     if (grid->globals->surface_type == LID_LOVE || grid->globals->surface_type == LID_MEMBR)
//     {
//      factor /= radius;
//      factor *= pow(radius + grid->globals->shell_thickness.Value(), 2.0)/radius;
//     }
//
//     sin2Lat = &(grid->trig2Lat(0,1));
//
//     for (i=0; i<node_num; i++) {
//         j = i*2;
//
//         soln(i,0) = 0.0;                    // no lon gradient
//         soln(i,1) = factor * sin2Lat[j];    // lat gradient
//     }
//
// };
//
//
// // -----------------------------------------------------------------------------
// // Full obliquity tidal forcing, degree-2 (e.g., Matsuyama 2014, Tyler 2011)
// // -----------------------------------------------------------------------------
// void deg2Obliq(Mesh * grid, Array2D<double> & soln, double simulationTime, double radius, double omega, double theta)
// {
//     double * sinLat, * sinLon, *cosLon, * cos2Lat;
//     double * sin2Lat;
//     int i,j, node_num;
//     double * val;
//     double * m;
//     double cosM, factor;
//
//     node_num = grid->node_num;
//
//     Array1D<double> * scalar_dummy;
//     scalar_dummy = new Array1D<double>(node_num);
//
//
//     // factor = pow(omega,2.0)*pow(radius,2.0)*ecc;
//     factor = -3. * grid->globals->loveReduct.Value() * pow(omega,2.0)*radius*theta;
//     cosM = cos(omega*simulationTime);
//     // sinM = sin(omega*simulationTime);
//
//
//     if (grid->globals->surface_type == LID_LOVE || grid->globals->surface_type == LID_MEMBR)
//     {
//      factor /= radius;
//      factor *= pow(radius + grid->globals->shell_thickness.Value(), 2.0)/radius;
//     }
//
//     // Assign pointers to start of trig node arrays
//     cos2Lat = &(grid->trig2Lat(0,0));
//     cosLon = &(grid->trigLon(0,0));
//     sinLon = &(grid->trigLon(0,1));
//     sinLat = &(grid->trigLat(0,1));
//     sin2Lat = &(grid->trig2Lat(0,1));
//
//     // Solve for dUdlon
//     for (i=0; i<node_num; i++) {
//         j = i*2;
//         // calculate potential gradient in longitude
//         // soln(i,0) += -factor * sinLat[j] * sinLon[j] * cosM;
//         //
//         //
//         // // calculate potential gradient in latitude
//         // soln(i,1) += factor * cos2Lat[j] * cosLon[j] * cosM;
//
//         // Calculate tidal potential
//         (*scalar_dummy)(i) = radius*factor/2. * cosM * sin2Lat[j] * cosLon[j];
//
//     }
//
//     pressureGradient(grid, soln, *scalar_dummy, node_num, -1.0);
//
//     delete scalar_dummy;
// };
//
//
// // -----------------------------------------------------------------------------
// // Westward travelling obliquity tidal forcing, degree-2 (e.g., Matsuyama 2014)
// // -----------------------------------------------------------------------------
// void deg2ObliqWest(Mesh * grid, Array2D<double> & soln, double simulationTime, double radius, double omega, double theta)
// {
//     double * sinLat, * sinLon, *cosLon, * cos2Lat;
//     int i,j, node_num;
//     double * val;
//     double * m;
//     double cosM, sinM, factor;
//
//     node_num = grid->node_num;
//
//     // factor = pow(omega,2.0)*pow(radius,2.0)*ecc;
//     factor = 1.5 * grid->globals->loveReduct.Value() * pow(omega,2.0)*radius*theta;
//     cosM = cos(omega*simulationTime);
//     sinM = sin(omega*simulationTime);
//
//     if (grid->globals->surface_type == LID_LOVE || grid->globals->surface_type == LID_MEMBR)
//     {
//      factor /= radius;
//      factor *= pow(radius + grid->globals->shell_thickness.Value(), 2.0)/radius;
//     }
//
//     // Assign pointers to start of trig node arrays
//     cos2Lat = &(grid->trig2Lat(0,0));
//     cosLon = &(grid->trigLon(0,0));
//     sinLon = &(grid->trigLon(0,1));
//     sinLat = &(grid->trigLat(0,1));
//
//     // Solve for dUdlon
//     for (i=0; i<node_num; i++) {
//         j = i*2;
//         // calculate potential gradient in longitude
//         soln(i,0) = -factor * sinLat[j] * (sinM * cosLon[j] + sinLon[j] * cosM);
//
//
//         // calculate potential gradient in latitude
//         soln(i,1) = factor * cos2Lat[j] * (cosM * cosLon[j] - sinM * sinLon[j]);
//
//     }
// };
//
//
// // -----------------------------------------------------------------------------
// // Eastward travelling obliquity tidal forcing, degree-2 (e.g., Matsuyama 2014)
// // -----------------------------------------------------------------------------
// void deg2ObliqEast(Mesh * grid, Array2D<double> & soln, double simulationTime, double radius, double omega, double theta)
// {
//     double * sinLat, * sinLon, *cosLon, * cos2Lat;
//     int i,j, node_num;
//     double * val;
//     double * m;
//     double cosM, sinM, factor;
//
//     node_num = grid->node_num;
//
//     // factor = pow(omega,2.0)*pow(radius,2.0)*ecc;
//     factor = 1.5 * grid->globals->loveReduct.Value() * pow(omega,2.0)*radius*theta;
//
//     if (grid->globals->surface_type == LID_LOVE || grid->globals->surface_type == LID_MEMBR)
//     {
//      factor /= radius;
//      factor *= pow(radius + grid->globals->shell_thickness.Value(), 2.0)/radius;
//     }
//
//     cosM = cos(omega*simulationTime);
//     sinM = sin(omega*simulationTime);
//
//     // Assign pointers to start of trig node arrays
//     cos2Lat = &(grid->trig2Lat(0,0));
//     cosLon = &(grid->trigLon(0,0));
//     sinLon = &(grid->trigLon(0,1));
//     sinLat = &(grid->trigLat(0,1));
//
//     for (i=0; i<node_num; i++) {
//         j = i*2;
//
//         // calculate potential gradient in longitude
//         soln(i,0) = -factor * sinLat[j] * (cosM * sinLon[j] - sinM * cosLon[j]);
//
//
//         // calculate potential gradient in latitude
//         soln(i,1) = factor * cos2Lat[j] * (sinM * sinLon[j] + cosM * cosLon[j]);
//
//     }
// };
//
//
// // -----------------------------------------------------------------------------
// // Full time-varying tidal forcing (ecc + obliq), degree-2 (e.g., Matsuyama 2014)
// // -----------------------------------------------------------------------------
// void deg2Full(Mesh * grid, Array2D<double> & soln, double simulationTime, double radius, double omega, double theta, double ecc)
// {
//     double cosM, sinM, factor_lon, factor_lat;              // cos(Mean anomaly), sin(Mean anomaly)
//     double * sinLon, * cosLon;
//     double * sinLat, * cosLat;
//     double * sin2Lat, * cos2Lat;
//     double * sin2Lon, * cos2Lon;
//
//     int i,j, node_num;
//     double * val;
//     double * m;
//
//     node_num = grid->node_num;
//
//     Array1D<double> * scalar_dummy;
//     scalar_dummy = new Array1D<double>(node_num);
//     // factor = pow(omega,2.0)*pow(radius,2.0)*ecc;
//
//     factor_lon = grid->globals->loveReduct.Value() * pow(omega,2.0)*radius;
//
//     if (grid->globals->surface_type == LID_LOVE || grid->globals->surface_type == LID_MEMBR)
//     {
//      factor_lon /= radius;
//      factor_lon *= pow(radius + grid->globals->shell_thickness.Value(), 2.0)/radius;
//     }
//
//     factor_lat = factor_lon;
//
//
//
//     cosM = cos(omega*simulationTime);
//     sinM = sin(omega*simulationTime);
//
//     // Assign pointers to start of trig node arrays
//     cosLon = &(grid->trigLon(0,0));
//     sinLon = &(grid->trigLon(0,1));
//     cosLat = &(grid->trigLat(0,0));
//     sinLat = &(grid->trigLat(0,1));
//
//     cos2Lat = &(grid->trig2Lat(0,0));
//     sin2Lat = &(grid->trig2Lat(0,1));
//     cos2Lon = &(grid->trig2Lon(0,0));
//     sin2Lon = &(grid->trig2Lon(0,1));
//
//     // Solve for dUdlon
//     // for (i=0; i<node_num; i++) {
//     //     j = i*2;
//     //     // calculate potential gradient in longitude
//     //     soln(i,0) = factor_lon * (3. * theta * sinLat[j] * sinLon[j] * cosM
//     //                     + 1.5 * ecc * cosLat[j]
//     //                     * (4. * sinM * cos2Lon[j]
//     //                     - 3. * cosM * sin2Lon[j]));
//     //
//     //     // calculate potential gradient in latitude
//     //     soln(i,1) = -factor_lat * (3. * theta * cos2Lat[j] * cosLon[j] * cosM
//     //                     + 0.75 * ecc * sin2Lat[j] * (3. * cosM * (1 + cos2Lon[j])
//     //                     + 4. * sin2Lon[j] * sinM));
//     // }
//
//
//     for (i=0; i<node_num; i++) {
//         j = i*2;
//         // calculate potential gradient in longitude
//         // soln(i,0) += -factor * sinLat[j] * sinLon[j] * cosM;
//         //
//         //
//         // // calculate potential gradient in latitude
//         // soln(i,1) += factor * cos2Lat[j] * cosLon[j] * cosM;
//
//         // Calculate obliquity tidal potential
//         (*scalar_dummy)(i) = -1.5*radius*theta*factor_lon * cosM * sin2Lat[j] * cosLon[j];
//
//         // Calculate eccentricity tidal potential
//         (*scalar_dummy)(i) += 0.75*factor_lon*ecc*radius *((1-3*sinLat[j]*sinLat[j])*cosM
//                                 + cosLat[j]*cosLat[j] * (3*cosM*cos2Lon[j]
//                                 + 4*sinM*sin2Lon[j]));
//
//     }
//
//     pressureGradient(grid, soln, *scalar_dummy, node_num, -1.0);
//
//     delete scalar_dummy;
// };
//
//
// // -----------------------------------------------------------------------------
// // Static + time-varying eccentricity tidal forcing, degree-2 (e.g., Matsuyama 2014)
// // -----------------------------------------------------------------------------
// void deg2EccTotal(Mesh * grid, Array2D<double> & soln, double simulationTime, double radius, double omega, double theta, double ecc)
// {
//     double cosM, sinM, factor_lon, factor_lat;              // cos(Mean anomaly), sin(Mean anomaly)
//     double * cosLat;
//     double * sin2Lat;
//     double * sin2Lon, * cos2Lon;
//     double factor;
//
//     int i,j, node_num;
//     double * val;
//     double * m;
//
//     node_num = grid->node_num;
//     //
//     // factor_lon = 0.5 * grid->globals->loveReduct.Value() * pow(omega,2.0)*radius;
//     //
//     // if (grid->globals->surface_type == LID_LOVE || grid->globals->surface_type == LID_MEMBR)
//     // {
//     //  factor_lon /= radius;
//     //  factor_lon *= pow(radius + grid->globals->shell_thickness.Value(), 2.0)/radius;
//     // }
//
//     cosM = cos(omega*simulationTime);
//     sinM = sin(omega*simulationTime);
//     //
//     // factor_lon *= (1. + 3.*ecc*cosM);
//     // factor_lat = factor_lon;
//     //
//     // // Assign pointers to start of trig node arrays
//     cosLat = &(grid->trigLat(0,0));
//
//     sin2Lat = &(grid->trig2Lat(0,1));
//     cos2Lon = &(grid->trig2Lon(0,0));
//     sin2Lon = &(grid->trig2Lon(0,1));
//
//     // Solve for dUdlon
//     // for (i=0; i<node_num; i++) {
//     //     j = i*2;
//     //     // factor_lon = factor_lat/cosLat[j];
//     //     // if (i < 2) factor_lon = 0.0;
//     //     // calculate potential gradient in longitude
//     //     soln(i,0) = factor_lon * (-3. * cosLat[j] * sin2Lon[j]
//     //                                   + 12. * ecc * cosLat[j] * cos2Lon[j] * sinM);
//     //
//     //     // calculate potential gradient in latitude
//     //     soln(i,1) = factor_lat * (-3. * sin2Lat[j] * cosSqLon[j]
//     //                                   - 12. * ecc * sin2Lat[j] * sinLon[j] * cosLon[j] * sinM);
//     // }
//
//     factor = grid->globals->loveReduct.Value() * pow(omega,2.0)*radius*ecc;
//
//     if (grid->globals->surface_type == LID_LOVE || grid->globals->surface_type == LID_MEMBR)
//     {
//      factor /= radius;
//      factor *= pow(radius + grid->globals->shell_thickness.Value(), 2.0)/radius;
//     }
//
//     for (i=0; i<node_num; i++) {
//         j = i*2;
//
//         soln(i,0) = factor * 1.5 * cosLat[j]   // cosLat here as cos^2(lat)/cos(lat)
//                         * (4.*sinM * cos2Lon[j]
//                         - 3.*cosM * sin2Lon[j]);
//
//         soln(i,0) += factor * 1.5 * cosLat[j]* sin2Lon[j];
//
//         soln(i,1) = -factor*0.75*sin2Lat[j]
//                         *(3.*cosM*(1.+cos2Lon[j])
//                         + 4.*sinM*sin2Lon[j]);
//
//         soln(i,1) += factor * 0.75 * sin2Lat[j] * ( 1. + cos2Lon[j]);
//
//     }
//
//     // avgAtPoles(grid, soln);
// };
//
// void deg2Planet(Mesh * grid, Array2D<double> & soln, Array1D<double> & potential, double simulationTime, double radius)
// {
//     double factor;              // cos(Mean anomaly), sin(Mean anomaly)
//     double * sinLon, * cosLon;
//     double * cosLat;
//
//     int i,j, node_num;
//     double * val;
//     double * m;
//
//     node_num = grid->node_num;
//
//     radius += grid->globals->shell_thickness.Value();
//
//     Array1D<double> * scalar_dummy;
//     scalar_dummy = new Array1D<double>(node_num);
//
//     double a1, a2;  //semimajor axes
//     double n1, n2;  //mean motions
//     double m2;
//     double cosnt, sinnt;
//     double cosphi, sinphi;
//     double dist;
//
//     // m2 = 8.931938e+22;  // Io
//     m2 = 1.4815e23;   // Ganymede
//     // m2 = 4.799e22;  //Europa
//
//
//     // a1 = 421800000.0; // Io
//     a2 = 1.074e9;   // Ganymede
//     a1 = grid->globals->a.Value();
//     // a2 = 671100000.0; // Europa
//
//     // n2 = 4.11e-5;   // Io
//     // n2 = 1.016e-5;     // Ganymede
//     // n2 = 2.05e-5;     // Europa
//     n1 = grid->globals->angVel.Value();
//     n2 = n1*0.5;
//
//     double nij = n1-n2;
//
//     cosnt = cos(nij*simulationTime);
//     sinnt = sin(nij*simulationTime);
//
//     // cosnt = cos(-n1*simulationTime);
//     // sinnt = sin(-n1*simulationTime);
//
//     dist = sqrt(pow(a1, 2.0) + pow(a2, 2.0) - 2.*a1*a2*cosnt);
//
//     cosphi = (a1 - a2*cosnt)/dist;
//     sinphi = a2*sinnt/dist;
//
//     factor = -1.0/grid->globals->g.Value() * 0.5 * 6.67408e-11 * m2 * pow(radius/dist, 2.0)/dist;
//
//     // Assign pointers to start of trig node arrays
//     cosLon = &(grid->trigLon(0,0));
//     sinLon = &(grid->trigLon(0,1));
//
//     cosLat = &(grid->trigLat(0,0));
//
//     double cosgam;
//     for (i=0; i<node_num; i++) {
//         j = i*2;
//         // Calculate obliquity tidal potential
//
//         cosgam = cosLat[j] * (cosLon[j]*cosphi + sinLon[j]*sinphi);
//         potential(i) = factor * (3. * cosgam*cosgam - 1.0);
//     }
//
//
//
//     // pressureGradient(grid, soln, *scalar_dummy, node_num, -1.0/grid->globals->g.Value());
//
//
//     // delete scalar_dummy;
// };
//
// void deg2General(Mesh * grid, Array2D<double> & soln, double simulationTime, double radius, Array1D<double> & scalar)
// {
//     double cos_qnt, sin_qnt;
//     double * cos2Lon, * sin2Lon;
//     double * cosLon, * cosLat, * sinLon;
//     double * P20, * P22;
//     double n1, n2, nij;
//     int i,j, node_num, q, q_max;
//
//
//     n1 = grid->globals->angVel.Value();
//     n2 = 0.5*n1;
//     nij = n1-n2;
//
//     radius += grid->globals->shell_thickness.Value();
//
//     node_num = grid->node_num;
//
//     cos2Lon = &(grid->trig2Lon(0,0));
//     sin2Lon = &(grid->trig2Lon(0,1));
//
//     cosLon = &(grid->trigLon(0,0));
//     sinLon = &(grid->trigLon(0,1));
//
//     cosLat = &(grid->trigLat(0,0));
//
//     P20 = &(grid->Pbar_20(0));
//     P22 = &(grid->Pbar_22(0));
//
//     Array1D<double> * scalar_dummy;
//     scalar_dummy = new Array1D<double>(node_num);
//
//     double a20, a22, b22;
//     q_max = grid->globals->freq.Value();
//
//     a20 = grid->globals->a20q[q_max];
//     a22 = grid->globals->a22q[q_max];
//     b22 = grid->globals->b22q[q_max];
//
//     // std::cout<<a20<<' '<<a22<<' '<<b22<<' '<<q<<std::endl;
//
//     cos_qnt = cos( (double)q_max * nij * simulationTime);
//     sin_qnt = sin( (double)q_max * nij * simulationTime);
//
//     double m2 = 4.799844e+22;  //Europa
//
//
//     // a2 = 421800000.0; // Io
//     // a2 = 1.074e9;   // Ganymede
//     double a1 = grid->globals->a.Value();
//     double a2 = 671100000.0; // Europa
//
//     // n2 = 4.11e-5;   // Io
//     // n2 = 1.016e-5;     // Ganymede
//     // n2 = 2.05e-5;     // Europa
//     // n1 = grid->globals->angVel.Value();
//     // n2 = n1*0.5;
//     //
//     // double nij = n1-n2;
//
//     double cosnt = cos(nij*(double)simulationTime);
//     double sinnt = sin(nij*(double)simulationTime);
//
//
//     // cosnt = cos(-n1*simulationTime);
//     // sinnt = sin(-n1*simulationTime);
//
//     double dist = sqrt(std::pow(a1, 2.0) + std::pow(a2, 2.0) - 2.*a1*a2*cosnt);
//
//     double cosphi = (a1 - a2*cosnt)/dist;
//     double sinphi = a2*sinnt/dist;
//
//     double factor = 0.5 * 6.67408e-11 * m2 * std::pow(radius/dist, 2.0)/dist;
//
//     // std::cout<<std::setprecision(8)<<std::scientific<<radius<<' '<<dist<<std::endl;
//     // Assign pointers to start of trig node arrays
//     cosLon = &(grid->trigLon(0,0));
//     sinLon = &(grid->trigLon(0,1));
//
//     cosLat = &(grid->trigLat(0,0));
//
//     double cosgam;
//     // for (i=0; i<node_num; i++) {
//     //     j = i*2;
//     //     // (*scalar_dummy)(i) = 0;
//     //
//     //     (*scalar_dummy)(i) = ((a20 * cos_qnt)* P20[i]
//     //                         + (a22 * cos_qnt) * P22[i] * cos2Lon[j]
//     //                         +b22 * P22[i] * sin2Lon[j] * sin_qnt)*grid->globals->loveReduct.Value();
//     // }
//
//
//     for (i=0; i<node_num; i++) {
//         j = i*2;
//         double U20=0.0;
//         double U22a=0.0;
//         double U22b=0.0;
//         for (q=q_max; q<=q_max; q++)
//         {
//             a20 = grid->globals->a20q[q];
//             a22 = grid->globals->a22q[q];
//             b22 = -grid->globals->b22q[q];
//
//             cos_qnt = cos( (double)q * nij * simulationTime);
//             sin_qnt = sin( (double)q * nij * simulationTime);
//
//             // (*scalar_dummy)(i) += ((a20 * cos_qnt)* P20[i]
//             //                         + (a22 * cos_qnt) * P22[i] * cos2Lon[j]
//             //                         +b22 * P22[i] * sin2Lon[j] * sin_qnt)*grid->globals->loveReduct.Value();
//
//
//             // std::cout<<q<<' '<<a20<<' '<<a22<<' '<<b22<<std::endl;
//             U20 += (a20 * cos_qnt);
//             U22a += (a22 * cos_qnt);
//             U22b += (b22 * sin_qnt);
//         }
//
//         U20 = (0.5*grid->globals->a20q[0] + U20)*P20[i];
//         U22a = (0.5*grid->globals->a22q[0] + U22a)*P22[i]*cos2Lon[j];
//         U22b = (-0.5*grid->globals->b22q[0] + U22b)*P22[i]*sin2Lon[j];
//
//         (*scalar_dummy)(i) = (U20 + U22a + U22b)*grid->globals->loveReduct.Value();
//
//         // cosgam = cosLat[j] * (cosLon[j]*cosphi + sinLon[j]*sinphi);
//         // double test = factor * (3. * pow(cosgam, 2.0) - 1.0)*grid->globals->loveReduct.Value();
//
//         // std::cout<<test<<' '<<(*scalar_dummy)(i)<<' '<<test/(*scalar_dummy)(i)<<std::endl;
//     }
//
//     // if (simulationTime <= 0.0+1e-5)
//     // {
//     //     for (i=0; i<node_num; i++) {
//     //         scalar(i) = grid->globals->loveReduct.Value()*(*scalar_dummy)(i)/(grid->globals->g.Value()*grid->globals->shell_factor_beta[2]);
//     //     }
//     // }
//
//     // grid->globals->Output->TerminateODIS();
//
//     pressureGradient(grid, soln, *scalar_dummy, node_num, -1.0);
//
//     delete scalar_dummy;
// };
//
// void deg2PlanetObl(Mesh * grid, Array2D<double> & soln, double simulationTime, double radius, double obl)
// {
//     double factor;              // cos(Mean anomaly), sin(Mean anomaly)
//     double * sinLon, * cosLon;
//     double * cosLat, * sinLat;
//
//     int i,j, node_num;
//     double * val;
//     double * m;
//
//     node_num = grid->node_num;
//
//
//     Array1D<double> * scalar_dummy;
//     scalar_dummy = new Array1D<double>(node_num);
//
//     double a1, a2;  //semimajor axes
//     double n1, n2;  //mean motions
//     double m2;
//     double cosnt, sinnt, cosM;
//     double cosphi, sinphi;
//     double dist;
//
//     // m2 = 8.931938e+22;  // Io
//     m2 = 1.4815e23;   // Ganymede
//     // m2 = 4.799e22;  //Europa
//
//
//     // a2 = 421800000.0; // Io
//     a1 = 1.074e9;   // Ganymede
//     a2 = 671100000.0; // Europa
//     // a1 = grid->globals->a.Value();
//     // a2 = 671100000.0; // Europa
//
//     // n2 = 4.11e-5;   // Io
//     // n2 = 1.016e-5;     // Ganymede
//     // n2 = 2.05e-5;     // Europa
//     n1 = grid->globals->angVel.Value();
//     n2 = 0.5*n1;
//
//     cosnt = cos((n1-n2)*simulationTime);
//     sinnt = sin((n1-n2)*simulationTime);
//     cosM  = cos(n1*simulationTime);
//
//     dist = sqrt(pow(a1, 2.0) + pow(a2, 2.0) - 2.*a1*a2*cosnt);
//
//     cosphi = (a1 - a2*cosnt)/dist;
//     sinphi = a2*sinnt/dist;
//
//     factor = 0.5 * 6.67e-11 * m2 * pow(radius/dist, 2.0)/dist;
//
//     // Assign pointers to start of trig node arrays
//     cosLon = &(grid->trigLon(0,0));
//     sinLon = &(grid->trigLon(0,1));
//
//     cosLat = &(grid->trigLat(0,0));
//     sinLat = &(grid->trigLat(0,1));
//
//     double cosgam;
//     for (i=0; i<node_num; i++) {
//         j = i*2;
//         // Calculate obliquity tidal potential
//         cosgam = cosLat[j] * sinLat[j] * (cosLon[j]*cosphi + sinLon[j]*sinphi) * cosphi*cosM;
//         (*scalar_dummy)(i) = factor * (3. * 2. * obl * cosgam - 1.0);
//     }
//
//     pressureGradient(grid, soln, *scalar_dummy, node_num, -1.0);
//
//     delete scalar_dummy;
// };
