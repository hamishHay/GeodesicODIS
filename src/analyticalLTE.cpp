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
#include "analyticalLTE.h"
#include "array2d.h"
#include <complex>
#include <array>
#include "mathRoutines.h"
#include "math.h"

const   std::complex<double> j(0.0,1.0); 

double pn(double n, double m) {
    return (n+1)*(n+m)/(n*(2*n+1));
}

double qn(double n, double m) {
    return n*(n+1-m)/((n+1)*(2*n+1));
}

std::complex<double> Kn(double n, double m, double lambs, double lam, double alpha, double Omega) {
    return lam + m/(n*(n+1)) - n*(n+1)/(lambs * lam) + j*alpha / (2*Omega);
}

std::complex<double> Ln(double n, double m, double lambs, double lam, double alpha, double Omega) {
    return lam + m/(n*(n+1)) + j*alpha / (2*Omega);
}


std::array<double, 6> analyticalLTE(Globals * consts, Mesh * grid, int forcing_type, double lat, double lon, double t)
{
    // int i,j, node_num, face_num;
    int i;
    // node_num = grid->node_num;
    // face_num = grid->face_num;

    // double factor, factor2;
    // double radius, omega;

    double * sinLon, * cosLon;
    double * sin2Lon, *cos2Lon;
    double * cosLat, * sinLat;
    double * sin2Lat, *cos2Lat;
    double * sinSqLat, * cosSqLat;

    // double cosM, sinM;
    // double cos2M, sin2M;

    double radius = consts->radius.Value();
    double Omega = consts->angVel.Value();
    double g = consts->g.Value();
    double h = consts->h.Value();
    double w = -Omega;
    double obl = consts->theta.Value();
    double lam = w/(2*Omega);
    double lambs = 4*Omega*Omega*radius*radius/(g*h);
    double alpha = consts->alpha.Value();
    double m = 1;



    // if (consts->surface_type == LID_LOVE || consts->surface_type == LID_MEMBR)
    // {
    //  radius += consts->shell_thickness.Value();
    // }


    // Assign pointers to start of trig node arrays
    cosLon = &(grid->trigLon(0,0));
    sinLon = &(grid->trigLon(0,1));
    cosLat = &(grid->trigLat(0,0));
    sinLat = &(grid->trigLat(0,1));


    // cos2Lon = &(grid->trig2Lon(0,0));
    // sin2Lon = &(grid->trig2Lon(0,1));
    // sin2Lat = &(grid->trig2Lat(0,1));

    // cos2Lat = &(grid->trig2Lat(0,0));

    cosSqLat = &(grid->trigSqLat(0, 0));
    // sinSqLat = &(grid->trigSqLat(0, 1));

    std::complex<double> U;
    std::complex<double> dUdt;
    std::complex<double> V;
    std::complex<double> dVdt;
    std::complex<double> ETA;
    std::complex<double> dETAdt;


    switch (forcing_type)
    {
    //   case ECC:
    //   {
    //     //   factor = 0.75 * consts->loveReduct.Value() * pow(omega,2.0)*pow(radius,2.0)*ecc;

    //     //   #pragma omp parallel for
    //     //   for (i=0; i<node_num; ++i) {
    //     //       j = i*2;

    //     //       potential(i) = factor*((1. - 3.*sinSqLat[i*2])*cosM
    //     //                               + cosSqLat[i*2] * (3.*cosM*cos2Lon[i*2]
    //     //                               + 4.*sinM*sin2Lon[i*2]));


    //     //     double U2 = 9./2. * cos2M * (0.75*cosSqLat[i*2] - 0.5);
    //     //     U2 += 17* 3/8. * cosSqLat[i*2] * (sin2M*sin2Lon[i*2] + cos2M*cos2Lon[2*i]);
    //     //     U2 *= factor * 4./3. * ecc;
    //     //     potential(i) += U2;

    //     //   }
    //   }
    //   break;

      case OBLIQ_WEST:
      {
        std::complex<double> L1 = Ln(1,1,lambs,lam,alpha,Omega);
        std::complex<double> K2 = Kn(2,1,lambs,lam,alpha,Omega);
        std::complex<double> U21 = 0.5*Omega*Omega * radius*radius * obl;
        std::complex<double> Psi11 = U21/(2*Omega) / (qn(1,1) - K2*L1/pn(2,1));
        std::complex<double> Phi21 = -j*L1/pn(2,1)*Psi11;
        // std::complex 

        std::complex<double> Y31 = 1.5 * (5*pow(cos(pi*0.5 - lat),2.0) - 1)*sin(pi*0.5 - lat)*exp(j * 1.0 * lon);
        std::complex<double> Y21 = 3 * cos(pi*0.5 - lat)*sin(pi*0.5 - lat) * std::exp(j * 1.0 * lon);
        std::complex<double> Y11 = sin(pi*0.5 - lat)*std::exp(j*1.0*lon);

        U = (-1/3. * (Psi11 * Y21) + j*Phi21*Y21)*std::exp(-j*w*t);
        dUdt = -j*w*U;
        U += std::conj(U);         // Add complex conjugate
        U /= 2*radius*sin(pi*0.5 - lat);

        dUdt += std::conj(dUdt);
        dUdt /= 2*radius*sin(pi*0.5 - lat);

        V = -(j*Psi11 *Y11 + Phi21*(4.0*Y31 - 9.0*Y11)/5.)*std::exp(-j*w*t);
        dVdt = -j*w*V;
        V += std::conj(V);         // Add complex conjugate 
        V /= 2*radius*sin(pi*0.5 - lat);

        // V *= -1.0;
        // dVdt *= -1.0;

        dVdt += std::conj(dVdt);
        dVdt /= 2*radius*sin(pi*0.5 - lat);

        ETA = j*6.0/pow(radius,2.0) * h/w * Phi21 * Y21  * std::exp(-j*w*t);
        dETAdt = -j*w*ETA;
        ETA += std::conj(ETA);
        ETA *= 0.5;

        dETAdt += std::conj(dETAdt);
        dETAdt *= 0.5;

        // double Uedge = grid->face_normal_vec_map(i, 0)*std::real(U) + grid->face_normal_vec_map(i, 1)*std::real(V);
        // double dUdtedge = grid->face_normal_vec_map(i, 0)*std::real(dUdt) + grid->face_normal_vec_map(i, 1)*std::real(dVdt);

        std::array<double, 6> data = {std::real(U), std::real(V), std::real(dUdt), std::real(dVdt), std::real(ETA), std::real(dETAdt)};

        return data;

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

