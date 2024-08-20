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

#include "mesh.h"
#include "globals.h"
#include "tidalPotentials.h"
#include "array2d.h"
#include "spatialOperators.h"
#include "gridConstants.h"
#include <math.h>
#include <iomanip>
#include <stdexcept>



void forcing(Globals * consts, Mesh * grid, Array1D<double> & potential, int forcing_type, double time, double ecc, double obl)
{
    int i,j;

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
          for (i=0; i<NODE_NUM; ++i) {
              j = i*2;

              potential(i) = factor*((1. - 3.*sinSqLat[i*2])*cosM
                                      + cosSqLat[i*2] * (3.*cosM*cos2Lon[i*2]
                                      + 4.*sinM*sin2Lon[i*2]));


            //double U2 = 9./2. * cos2M * (0.75*cosSqLat[i*2] - 0.5);
            //U2 += 17* 3/8. * cosSqLat[i*2] * (sin2M*sin2Lon[i*2] + cos2M*cos2Lon[2*i]);
            //U2 *= factor * 4./3. * ecc;
            //potential(i) += U2;

          }
      }
      break;

      case OBLIQ:
      {
          factor = -3./2. * consts->loveReduct.Value() * pow(omega,2.0)*pow(radius,2.0)*obl;

          #pragma omp parallel for
          for (i=0; i<NODE_NUM; ++i) {
              j = i*2;

              potential(i) = factor * cosM * sin2Lat[j] * cosLon[j];
          }

      }
      break;

      case OBLIQ_WEST:
      {
          factor = 0.5 * consts->loveReduct.Value() * pow(omega,2.0)*pow(radius,2.0)*obl;

          #pragma omp parallel for
          for (i=0; i<NODE_NUM; ++i) {
              j = i*2;

              potential(i) = 3*factor * sinLat[j]*cosLat[j]*(cosLon[j]*cosM - sinLon[j]*sinM);
              
          }

      }
      break;

      case FULL2:
      {
          factor = 1/32. * consts->loveReduct.Value() * pow(omega,2.0)*pow(radius,2.0);
          for (i=0; i<NODE_NUM; i++) {
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
          for (i=0; i<NODE_NUM; ++i) {
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
          for (i=0; i<NODE_NUM; ++i) {
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
            for (i=0; i<NODE_NUM; i++) {
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

        case GIANT_IMPACT:
        {
            double M_impactor = 6.4171E23;
            double omegaT = omega*time;
            double cosphi, sinphi;
            double rx, ry, rmag;
            double t = time/86400.0;
            double cosGam;
            double v_c = 1.0;       // Velocity at contact, in units of escape velocity 
            double b = 0.7;         // Impact parameter
            double t_to_impact = 40*3600.0 - time;

            // t_to_impact = 40*3600.0 - time;
            try {
                // impact_pos_vel_b_v_c_t(rx, ry, t_to_impact, b, v_c, consts->radius.Value(), 3000e3, 6e24, M_impactor);
                impact_pos_vel_b_v_c_t(rx, ry, t_to_impact, b, v_c, consts->radius.Value(), 3396.2e3, 5.972e24, M_impactor);
            }
            catch (const std::exception& ex) {
                // impactor has impacted!
                consts->Output->TerminateODIS();
            }
            
            rmag = sqrt(rx*rx + ry*ry);


            cosphi = (rx*cosM + ry*sinM)/rmag;
            sinphi = (-rx*sinM + ry*cosM)/rmag;


            double fac = 6.67e-11*M_impactor/rmag * pow(radius/rmag, 2.0); 
            for (i=0; i<NODE_NUM; i++) {
                j = i*2;

                cosGam = cosLat[j]*(cosLon[j]*cosphi + sinLon[j]*sinphi);

                potential(i) = fac * (3*cosGam*cosGam - 1.0)*0.5
                                + fac * (radius/rmag) * (5*pow(cosGam,3.0) - 3*cosGam)*0.5;

            }
            
        }
        break;

        case NONE:
            break;
    }
};


// Function ported from WOMA https://pypi.org/project/woma/
void impact_pos_vel_b_v_c_r(double &x_, double &y_, double &t_out_, double b, double v_c, double r, double R_t, double R_i, double M_t, double M_i)
{

    const double G = 6.67430e-11; // Gravitational constant

    double mu = G * (M_t + M_i);
    double v_esc = sqrt(2 * mu / (R_t + R_i));
    double r_c = R_t + R_i;

    // Convert to b and v_c (m/s) if necessary
    v_c *= v_esc;
   

    // Contact position
    double y_c = b * r_c;
    if (r < r_c) {
        throw std::invalid_argument(
            "Invalid r = " + std::to_string(r) + " m for body radii " + 
            std::to_string(R_t) + " + " + std::to_string(R_i) + " = " + std::to_string(r_c) + " m"
        );
    }

    double v, x, y, theta_c, theta, alpha_;
    double a, e;
    if (v_c == v_esc) { // Parabola
        v = sqrt(2 * mu / r);
        y = v_c * y_c / v;
        x = sqrt(r * r - y * y);

        theta_c = pi - std::acos(y_c * y_c * v_c * v_c / (mu * r_c) - 1);
        theta = pi - std::acos(y_c * y_c * v_c * v_c / (mu * r) - 1);

        alpha_ = theta_c / 2.0;
    } else { // Ellipse or hyperbola
        a = 1.0 / (2.0 / r_c - v_c * v_c / mu);
        v = sqrt(mu * (2.0 / r - 1.0 / a));
        y = v_c * y_c / v;
        x = sqrt(r * r - y * y);

        double r_p = std::min(
            fabs(a + sqrt(a * a - a * v_c * v_c * y_c * y_c / mu)),
            fabs(a - sqrt(a * a - a * v_c * v_c * y_c * y_c / mu))
        );
        e = 1.0 - r_p / a;

        double r_a = 2.0 * a - r_p;
        if (v_c < v_esc && r > r_a) {
            throw std::runtime_error(
                "Invalid r = " + std::to_string(r) + " m for bound orbit (v_esc = " + 
                std::to_string(v_esc) + " m/s) with apoapsis = " + std::to_string(r_a) + " m"
            );
        }

        theta_c = acos((1.0 - a * (1.0 - e * e) / r_c) / e);
        theta = acos((1.0 - a * (1.0 - e * e) / r) / e);

        alpha_ = asin(sqrt(a * a * (1.0 - e * e) / (2.0 * a * r_c - r_c * r_c)));
    }


    double t, t_c;
    if (b == 0.0) {
        if (v_c == v_esc) { // Parabolic
            t = sqrt(2.0 * r * r * r / (9.0 * mu));
            t_c = sqrt(2.0 * r_c * r_c * r_c / (9.0 * mu));
        } else if (a > 0) { // Elliptical
            double w = 1.0 / r - v * v / (2.0 * mu);
            double w_c = 1.0 / r_c - v_c * v_c / (2.0 * mu);
            double wr = w * r;
            double wr_c = w_c * r_c;
            t = (asin(sqrt(wr)) - sqrt(wr * (1.0 - wr))) / sqrt(2.0 * mu * w * w * w);
            t_c = (asin(sqrt(wr_c)) - sqrt(wr_c * (1.0 - wr_c))) / sqrt(2.0 * mu * w_c * w_c * w_c);
        } else { // Hyperbolic
            double w = abs(1.0 / r - v * v / (2.0 * mu));
            double w_c = abs(1.0 / r_c - v_c * v_c / (2.0 * mu));
            double wr = w * r;
            double wr_c = w_c * r_c;
            t = (sqrt(wr * wr + wr) - log(sqrt(wr) + sqrt(1.0 + wr))) / (sqrt(2.0 * mu * w * w * w));
            t_c = (sqrt(wr_c * wr_c + wr_c) - log(sqrt(wr_c) + sqrt(1.0 + wr_c))) / (sqrt(2.0 * mu * w_c * w_c * w_c));
        }
    } else {
        if (v_c == v_esc) { // Parabolic
            double E = tan(0.5 * (pi - theta));
            double E_c = tan(0.5 * (pi - theta_c));
            double M = E + E * E * E / 3.0;
            double M_c = E_c + E_c * E_c * E_c / 3.0;
            double r_p = mu * (1.0 + cos(pi - theta)) / (v * v);
            t = sqrt(2.0 * r_p * r_p * r_p / mu) * M;
            t_c = sqrt(2.0 * r_p * r_p * r_p / mu) * M_c;
        } else if (a > 0) { // Elliptical
            double E = acos((e + cos(pi - theta)) / (1.0 + e * cos(pi - theta)));
            double E_c = acos((e + cos(pi - theta_c)) / (1.0 + e * cos(pi - theta_c)));
            double M = E - e * sin(E);
            double M_c = E_c - e * sin(E_c);
            t = sqrt(a * a * a / mu) * M;
            t_c = sqrt(a * a * a / mu) * M_c;
        } else { // Hyperbolic
            double E = acosh((e + cos(pi - theta)) / (1.0 + e * cos(pi - theta)));
            double E_c = acosh((e + cos(pi - theta_c)) / (1.0 + e * cos(pi - theta_c)));
            double M = -E + e * sinh(E);
            double M_c = -E_c + e * sinh(E_c);
            t = sqrt(-(a * a * a) / mu) * M;
            t_c = sqrt(-(a * a * a) / mu) * M_c;
        }
    }
    t_out_ = t - t_c;

    double phi;
    if (b == 0.0) {
        phi = 0.0;
    } else {
        phi = alpha_ - theta_c + theta - asin(y / r);
    }

    x_ = x * cos(phi) - y * sin(phi);
    y_ = x * sin(phi) + y * cos(phi);
};

// Function ported from WOMA https://pypi.org/project/woma/
void impact_pos_vel_b_v_c_t(double &x_, double &y_, double t, double b, double v_c, double R_t, double R_i, double M_t, double M_i, double r_max_factor)
{
    double r_min = R_t + R_i;
    double r_max = r_min * r_max_factor;
    double r;

    int i = 0;
    const int i_max = 100;
    double t_ = 0.0;
    const double tol = 1e-6;
    double t_out, x, y;
    
    // Bisection to find the separation to give the desired time to impact
    while (tol < fabs(t_ - t) / t) {
        r = 0.5 * (r_min + r_max);
    
        i++;

        try {
            
            impact_pos_vel_b_v_c_r(x, y, t_, b, v_c, r, R_t, R_i, M_t, M_i);
        } catch (const std::exception& ex) {
            t_ = t * 2;  // Set t_ to a value ensuring r_max is reduced
            if (i >= i_max) {
                throw std::runtime_error("Failed to find r(t) after " + std::to_string(i) + " iterations");
            }
        }

        // Bisect
        if (t_ < t) {
            r_min = r;
        } else {
            r_max = r;
        }

        if (i >= i_max) {
            throw std::runtime_error("Failed to find r(t) after " + std::to_string(i) + " iterations");
        }
    }

    // update x_ and y_ initial position and velocity
    impact_pos_vel_b_v_c_r(x_, y_, t_out, b, v_c, r, R_t, R_i, M_t, M_i);
};

