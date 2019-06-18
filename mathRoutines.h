#ifndef MATHROUTINES_H_INCDLUDED
#define MATHROUTINES_H_INCDLUDED

// #include "field.h"
// #include "solver.h"
//#include "mesh.h"
#include <math.h>
#include <iostream>

const double pi = 3.1415926535897932384626433832795028841971693993751058;
const double radConv = pi / 180.0; // global constant for converting deg --> rad

inline void sph2cart(double &, double &, double &, double, double, double);
inline void sph2cart(double &x, double &y, double &z, double r, double colat, double lon)
{
    x = r*sin(colat)*cos(lon);
    y = r*sin(colat)*sin(lon);
    z = r*cos(colat);
}

// Mapping function from Lee and Macdonald (2009) to get to GSTC
// Returns a mapping factor that is used appy curvature to the local system of
// equations.
inline void mapFactorAtPoint(double &, double, double, double, double);
inline void mapFactorAtPoint(double &m, double lat1, double lat2, double lon1, double lon2)
{
  m = 2.0 / (1.0 + sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(lon2 - lon1));
}

// Mapping function from Lee and Macdonald (2009) to get to GSTC.
// Function assigns a mapping factor that is used appy curvature to the local system of
// equations, as well as the x and y coordinates in the local GSTC coordinate system.
// lat1 and lat2 are the coordinates at which the projection is centered around.
inline void mapAtPoint(double &, double &, double &, double &, double &, double &, double &, double &);
inline void mapAtPoint(double &m, double &x, double &y, double &lat1, double &lat2, double &lon1, double &lon2, double &r)
{
  mapFactorAtPoint(m, lat1, lat2, lon1, lon2);

  x = m * r * (cos(lat2)*sin(lon2-lon1));

  y = m * r * (sin(lat2)*cos(lat1) - cos(lat2)*sin(lat1)*cos(lon2-lon1));
}

// Velocity transform factors from Lee and Macdonald (2009) to get to GSTC.
// Function calculates the cos(alpha) and sin(alpha) transform coefficients,
// which are required to convert spherical coord velocities to mapped velocities
// and back again. The back conversion requires inversion of the tranform matrix:
//
//                            ( cos_a   sin_a)
//                        T = (-sin_a   cos_a)
//
inline void velTransform(double &, double &, double &, double &, double &, double &);
inline void velTransform(double &cos_a, double &sin_a, double &lat1, double &lat2, double &lon1, double &lon2)
{
  cos_a = cos(lat1) * cos(lat2);
  cos_a += (1.0 + sin(lat1)*sin(lat2))*cos(lon2-lon1);
  cos_a /= (1.0 + sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(lon2-lon1));

  sin_a = -(sin(lat1) + sin(lat2))*sin(lon2-lon1);
  sin_a /= (1.0 + sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(lon2-lon1));
}

inline void volumeSphericalTriangle(double &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double);
inline void volumeSphericalTriangle(double &vol, double &x1, double &x2, double &x3, double &y1, double &y2, double &y3, double &z1, double &z2, double &z3, double r)
{
    double bXc_x, bXc_y, bXc_z;
    double vTop;
    double magA, magB, magC;
    double AdotB, CdotA, BdotC;
    double solidAngle;

    bXc_x = y2*z3 - z2*y3;
    bXc_y = -(x2*z3 - z2*x3);
    bXc_z = x2*y3 - x3*y2;

    vTop = x1*bXc_x + y1*bXc_y + z1*bXc_z ;

    magA = sqrt(x1*x1 + y1*y1 + z1*z1);
    magB = sqrt(x2*x2 + y2*y2 + z2*z2);
    magC = sqrt(x3*x3 + y3*y3 + z3*z3);

    AdotB = x1*x2 + y1*y2 + z1*z2;
    CdotA = x3*x1 + y3*y1 + z3*z1;
    BdotC = x2*x3 + y2*y3 + z2*z3;

    solidAngle = 2.*atan2(vTop,(magA*magB*magC + AdotB*magC + CdotA*magB + BdotC*magA));


    vol = solidAngle/3.0 * pow(r, 3.0);
}

inline void distanceBetween(double &, double &, double &, double &, double &);
inline void distanceBetween(double &len, double &x1, double &x2, double &y1, double &y2)
{
  double xx, yy;

  xx = x2 - x1;
  yy = y2 - y1;

  len = sqrt(xx*xx + yy*yy);
}

inline void distanceBetweenSph(double &, double &, double &, double &, double &, double &);
inline void distanceBetweenSph(double &len, double &lat1, double &lat2, double &lon1, double &lon2, double &r)
{
  len = r * acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon2-lon1));
}

inline void midpointBetween(double &, double &, double &, double &, double &, double &);
inline void midpointBetween(double &xc, double &yc, double &x1, double &x2, double &y1, double &y2)
{
  xc = 0.5 * (x1 + x2);
  yc = 0.5 * (y1 + y2);
}

inline void midpointBetweenSph(double sph1[], double sph2[], double sphc[]);
inline void midpointBetweenSph(double sph1[], double sph2[], double sphc[])
{
  double xyz1[3], xyz2[3];
  double xyzc[3];

  xyz1[0] = cos(sph1[0])*cos(sph1[1]);
  xyz1[1] = cos(sph1[0])*sin(sph1[1]);
  xyz1[2] = sin(sph1[0]);

  xyz2[0] = cos(sph2[0])*cos(sph2[1]);
  xyz2[1] = cos(sph2[0])*sin(sph2[1]);
  xyz2[2] = sin(sph2[0]);

  for (int k=0; k<3; k++) {
      xyzc[k] = 0.5*(xyz1[k]+xyz2[k]);
  }

  sphc[0] = atan2(xyzc[2], sqrt(pow(xyzc[0], 2.0) + pow(xyzc[1], 2.0)));
  sphc[1] = atan2(xyzc[1], xyzc[0]);
}

// Function to find the outward unit vector between two points with coordinates
// p1 = (x1,y1) and p2 = (x2,y2)
inline void normalVectorBetween(double &, double &, double &, double &, double &, double &);
inline void normalVectorBetween(double &xn, double &yn, double &x1, double &x2, double &y1, double &y2)
{
  double xx,yy,mag;

  xx = x2-x1;                   // Components of the tangent vector
  yy = y2-y1;

  mag = sqrt(xx*xx + yy*yy);    // Find magnitude of tangent vector

  xx /= mag;                    // Normalise tangent vector
  yy /= mag;

  xn = -yy;                     // Find normal unit vector to the tangent vector
  yn = xx;
}

inline void normalVectorBetweenMap(double n[], double sph1[], double sph2[]);
inline void normalVectorBetweenMap(double n[], double sph1[], double sph2[])
{
  double xx,yy, c1[2], c2[2], mag, sph_mid[3];

  midpointBetweenSph(sph1, sph2, sph_mid);

  double latc = sph_mid[0];
  double lonc = sph_mid[1];
  double lat1 = sph1[0];
  double lon1 = sph1[1];
  double m = 0.0;
  double r = 1.0;

  mapFactorAtPoint(m, latc, lat1, lonc, lon1);
  c1[0] = m * r * (cos(lat1)*sin(lon1-lonc));
  c1[1] = m * r * (sin(lat1)*cos(latc) - cos(lat1)*sin(latc)*cos(lon1-lonc));

  lat1 = sph2[0];
  lon1 = sph2[1];

  mapFactorAtPoint(m, latc, lat1, lonc, lon1);
  c2[0] = m * r * (cos(lat1)*sin(lon1-lonc));
  c2[1] = m * r * (sin(lat1)*cos(latc) - cos(lat1)*sin(latc)*cos(lon1-lonc));

  xx = c2[0]-c1[0];                   // Components of the tangent vector
  yy = c2[1]-c1[1];

  mag = sqrt(xx*xx + yy*yy);    // Find magnitude of tangent vector

  xx /= mag;                    // Normalise tangent vector
  yy /= mag;

  n[0] = -yy;                     // Find normal unit vector to the tangent vector
  n[1] = xx;
}

// Function to find the outward unit vector between two points with coordinates
// p1 = (x1,y1) and p2 = (x2,y2)
inline void normalVectorBetweenSph(double &, double &, double &);
inline void normalVectorBetweenSph(double sph1[], double sph2[], double sph3[])
{
    double c1[3], c2[3], n[3], sph_mid[3];



    // Get cartesian coords of each spherical coordiate
    sph2cart(c1[0], c1[1], c1[2], 1.0, sph1[0], sph1[1]);
    sph2cart(c2[0], c2[1], c2[2], 1.0, sph2[0], sph2[1]);

    // Take cross product to find vector tangent to line between sph1 and sph2
    n[0] = c1[1]*c2[2] - c1[2]*c2[1];
    n[1] = -(c1[0]*c2[2] - c1[2]*c2[0]);
    n[2] = c1[0]*c2[1] - c1[1]*c2[0];

    double mag = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);

    n[0] /= mag;
    n[1] /= mag;
    n[2] /= mag;

    midpointBetweenSph(sph1, sph2, sph_mid);

    double lat = sph_mid[0];
    double lon = sph_mid[1];

    sph3[0] = sin(lat)*cos(lon)*n[0] + sin(lat)*sin(lon)*n[1] - cos(lat)*n[2];
    sph3[1] = -sin(lon)*n[0] + cos(lon)*n[1];
    sph3[2] = cos(lat)*cos(lon)*n[0] + cos(lat)*sin(lon)*n[1] + sin(lat)*n[2];

}

inline void intersectPointSph(double &, double &, double &, double &, double &);
inline void intersectPointSph(double sph1[], double sph2[], double sph3[], double sph4[], double sph_int[])
{
    double c1[3], c2[3], c3[3], c4[4], n1[3], n2[3], mp[3];

    // sph1[0] = pi*0.5 - sph1[0];
    // sph2[0] = pi*0.5 - sph2[0];
    // sph3[0] = pi*0.5 - sph3[0];
    // sph4[0] = pi*0.5 - sph4[0];
    //
    // // Get cartesian coords of each spherical coordiate
    // sph2cart(c1[0], c1[1], c1[2], 1.0, sph1[0], sph1[1]);
    // sph2cart(c2[0], c2[1], c2[2], 1.0, sph2[0], sph2[1]);
    // sph2cart(c3[0], c3[1], c3[2], 1.0, sph3[0], sph3[1]);
    // sph2cart(c4[0], c4[1], c4[2], 1.0, sph4[0], sph4[1]);

    c1[0] = cos(sph1[0])*cos(sph1[1]);
    c2[0] = cos(sph2[0])*cos(sph2[1]);
    c3[0] = cos(sph3[0])*cos(sph3[1]);
    c4[0] = cos(sph4[0])*cos(sph4[1]);

    c1[1] = cos(sph1[0])*sin(sph1[1]);
    c2[1] = cos(sph2[0])*sin(sph2[1]);
    c3[1] = cos(sph3[0])*sin(sph3[1]);
    c4[1] = cos(sph4[0])*sin(sph4[1]);

    c1[2] = sin(sph1[0]);
    c2[2] = sin(sph2[0]);
    c3[2] = sin(sph3[0]);
    c4[2] = sin(sph4[0]);



    // Take cross product to find vector tangent to line between sph1 and sph2
    n1[0] = c1[1]*c2[2] - c1[2]*c2[1];
    n1[1] = -(c1[0]*c2[2] - c1[2]*c2[0]);
    n1[2] = c1[0]*c2[1] - c1[1]*c2[0];

    n2[0] = c3[1]*c4[2] - c3[2]*c4[1];
    n2[1] = -(c3[0]*c4[2] - c3[2]*c4[0]);
    n2[2] = c3[0]*c4[1] - c3[1]*c4[0];

    double mag1 = sqrt(n1[0]*n1[0] + n1[1]*n1[1] + n1[2]*n1[2]);
    double mag2 = sqrt(n2[0]*n2[0] + n2[1]*n2[1] + n2[2]*n2[2]);

    n1[0] /= mag1;
    n1[1] /= mag1;
    n1[2] /= mag1;

    n2[0] /= mag2;
    n2[1] /= mag2;
    n2[2] /= mag2;

    mp[0] = n1[1]*n2[2] - n2[1]*n1[2];
    mp[1] = n2[0]*n1[2] - n1[0]*n2[2];
    mp[2] = n1[0]*n2[1] - n2[0]*n1[1];

    mag1 = sqrt(mp[0]*mp[0] + mp[1]*mp[1] + mp[2]*mp[2]);
    mp[0] /= mag1;
    mp[1] /= mag1;
    mp[2] /= mag1;

    sph_int[0] = asin(mp[2]/1.0);
    sph_int[1] = atan2(mp[1], mp[0]);// + pi;
    sph_int[2] = 1.0;

    if (sph_int[1] < 0.0) sph_int[1] += 2*pi;
}

inline void angleBetweenMap(double sph1[], double sph2[], double sphc[], double &angle);
inline void angleBetweenMap(double sph1[], double sph2[], double sphc[], double &angle)
{
    double latc = sphc[0];
    double lonc = sphc[1];
    double lat1 = sph1[0];
    double lon1 = sph1[1];
    double lat2 = sph2[0];
    double lon2 = sph2[1];
    double m = 0.0, x1=0.0, y1=0.0, x2=0.0, y2=0.0;
    double xc=0.0, yc=0.0;
    double r = 1.0;

    mapAtPoint(m, x1, y1, latc, lat1, lonc, lon1, r);
    mapAtPoint(m, x2, y2, latc, lat2, lonc, lon2, r);

    double dot = x1*x2 + y1*y2;     // dot product between [x1, y1] and [x2, y2]
    double det = x1*y2 - y1*x2;     // determinant
    angle = atan2(det, dot);

    double ang1 = atan2(y1 - yc, x1 - xc);
    double ang2 = atan2(y2 - yc, x2 - xc);
    angle = ang1 - ang2;

    if (angle > pi) angle -= 2*pi;
    else if (angle < -pi) angle += 2*pi;
}


// Function to find the area of a triangle with points p(xc,yc), p(x1,y2) and
// p(x2,y2) in cartesian or mapping coordinates.
inline void triangularArea(double &, double &, double &, double &, double &, double &, double &);
inline void triangularArea(double &area, double &xc, double &yc, double &x1, double &x2, double &y1, double &y2)
{
  // area = 0.5 * fabs((xc - x2)*(y1 - yc) - (xc - x1)*(y2 - yc));
  area = 0.5 * fabs(xc*(y1 - y2) + x1*(y2-yc) + x2*(yc-y1));
}

inline void triangularAreaSph(double &, double &, double &, double &, double &, double &, double &, double &);
inline void triangularAreaSph(double &area, double &lat1, double &lat2, double &lat3, double &lon1, double &lon2, double &lon3, double &r)
{

  double a, b, c;
  double A, B, C;

  // double a, b;

  // std::cout<<lat1<<' '<<lat2<<' '<<lat3<<std::endl;

  c = fabs(acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(fabs(lon2-lon1))));
  a = fabs(acos(sin(lat2) * sin(lat3) + cos(lat2) * cos(lat3) * cos(fabs(lon3-lon2))));
  b = fabs(acos(sin(lat3) * sin(lat1) + cos(lat3) * cos(lat1) * cos(fabs(lon1-lon3))));

  A = fabs(acos((cos(a) - cos(b)*cos(c))/(sin(b)*sin(c))));
  B = fabs(acos((cos(b) - cos(a)*cos(c))/(sin(a)*sin(c))));
  C = fabs(acos((cos(c) - cos(b)*cos(a))/(sin(b)*sin(a))));

  area = pow(r,2.0) * fabs((A + B + C) - pi);
}

//
// double arcLength(double lat1, double lat2, double lon1, double lon2);
//
// inline double gamm1n(double xx);
// inline double gamm1n(double xx) {
//     double x,y,tmp,ser;
//     static double cof[6]={76.18009172947146,
//                           -86.50532032941677,
//                           24.01409824083091,
//                           -1.231739572450155,
//                           0.1208650973866179e-2,
//                           -0.5395239384953e-5};
//     int j;
//
//     y=x=xx;
//     tmp=x+5.5;
//     tmp -= (x+0.5)*log(tmp);
//     ser=1.000000000190015;
//     for (j=0; j<=5; j++) ser += cof[j]/++y;
//     return -tmp+log(2.5066282746310005*ser/x);
// }
//
// inline double factrl(int n)
// {
//     double gammln(float xx);
//
//     static int ntop=4;
//     static float a[33]={1.0,
//                         1.0,
//                         2.0,
//                         6.0,
//                         24.0,
//                         120.0,
//                         720.0,
//                         5040.0,
//                         40320.0,
//                         362880.0,
//                         3628800.0,
//                         39916800.0,
//                         479001600.0,
//                         6227020800.0,
//                         87178291200.0,
//                         1307674368000.0,
//                         20922789888000.0,
//                         355687428096000.0,
//                         6402373705728000.0,
//                         1.21645100408832e+17,
//                         2.43290200817664e+18,
//                         5.109094217170944e+19,
//                         1.1240007277776077e+21,
//                         2.585201673888498e+22,
//                         6.204484017332394e+23,
//                         1.5511210043330984e+25,
//                         4.032914611266057e+26,
//                         1.0888869450418352e+28,
//                         3.048883446117138e+29,
//                         8.841761993739701e+30,
//                         2.652528598121911e+32,
//                         8.222838654177924e+33,
//                         2.6313083693369355e+35};
//
//     int j;
//
//     if (n > 32) return exp(gamm1n(n+1.0));
//     // while (ntop<n) {
//     //   j=ntop++;
//     //   a[ntop]=a[j]*ntop;
//     // }
//     return a[n];
// }
//
// inline double legendreNorm(int l, int m);
// inline double legendreNorm(int l, int m) {
//     double norm[9][9] = {{1.0, 0, 0, 0, 0, 0, 0, 0, 0},
//                          {1.73205080757, 1.73205080757, 0, 0, 0, 0, 0, 0, 0},
//                          {2.2360679775, 1.29099444874, 0.645497224368, 0, 0, 0, 0, 0, 0},
//                          {2.64575131106, 1.08012344973, 0.341565025532, 0.139443337756, 0, 0, 0, 0, 0},
//                          {3.0, 0.948683298051, 0.22360679775, 0.0597614304667, 0.0211288563682, 0, 0, 0, 0},
//                          {3.31662479036, 0.856348838578, 0.161834718743, 0.0330343736322, 0.00778627653585, 0.00246223683452, 0, 0, 0},
//                          {3.60555127546, 0.786795792469, 0.124403337882, 0.020733889647, 0.0037854730215, 0.000807065559928, 0.000232979759139, 0, 0},
//                          {3.87298334621, 0.731925054711, 0.0996023841112, 0.0140859042455, 0.00212352996432, 0.000353921660721, 6.94097482422e-05, 1.8550535516e-05, 0},
//                          {4.12310562562, 0.687184270936, 0.0821342300508, 0.0101100248374, 0.00130519859416, 0.000180998479074, 2.79286716592e-05, 5.09905448963e-06, 1.27476362241e-06}};
//
//     return norm[l][m];
// }
//
// // Numerical recipes function for calculating the associated Legendre
// // function of x.
// inline double assLegendre(int l, int m, double x) __attribute__((always_inline));
// inline double assLegendre(int l, int m, double x) {
//     double fact, pll, pmm, pmmp1, somx2, normalise;
//     double delta = 0.0;
//     int i, ll;
//
//     if (m==0) delta = 1.0;
//
//     // normalise = legendreNorm(l,m);
//     normalise =  sqrt((2.0-delta)*(2.0*(double)l+1.0)*factrl(l-m)/factrl(l+m));
//
//     // if (m%2 != 0) normalise = -normalise;
//
//     // std::cout<<"l="<<l<<" m="<<m<<" norm="<<normalise<<std::endl;
//
//     pmm = 1.0;
//     if (m > 0)
//     {
//         somx2 = sqrt((1.0 - x)*(1.0+x));
//         fact = 1.0;
//         for (i=1; i<=m; i++)
//         {
//             pmm *= -fact*somx2;
//             fact += 2.0;
//         }
//     }
//     if (l == m) {
//         return normalise*pmm;
//     }
//     else {
//         pmmp1 = x*(2*m+1)*pmm;
//         if (l == (m+1)) {
//             return normalise*pmmp1;
//         }
//         else {
//             for (ll=m+2; ll<=l; ll++) {
//                 pll = (x*(2*ll-1)*pmmp1 - (ll+m-1)*pmm)/(ll-m);
//                 pmm = pmmp1;
//                 pmmp1 = pll;
//             }
//             return normalise*pll;
//         }
//     }
// };
//
// inline double fastSin(double x);
// inline double fastSin(double x) {
//     double ans;
//     double xxx = x*x*x;
//     double xxxxx = xxx*x*x;
//
//     ans = x;
//     ans -= (xxx)/factrl(3);
//     ans += (xxxxx)/factrl(5);
//
//     return x;
// }
//
// //func returns value n at given i and j based on Lagrange Interpolation
// //of points n+1, n+2, and n+3 for a 2nd degree polynomial
// inline double linearInterp1Array(Field * field, double ** solution, int i, int j) __attribute__((always_inline));
// inline double linearInterp1Array(Field * field, double ** solution, int i, int j) {
//     //first order accurate forward difference method
//     if (i < field->fieldLatLen / 2) {
//         return solution[i + 1][j] + (solution[i + 1][j] - solution[i + 2][j]);
//     }
//     else {
//         return solution[i - 1][j] + (solution[i - 1][j] - solution[i - 2][j]);
//     }
// };
//
//
//
//
//
// double linearInterp2(Field * field, int i, int j);
//
// inline double lagrangeInterp2Array(Field * fieldDat, double ** field, int i, int j) __attribute__((always_inline));
// inline double lagrangeInterp2Array(Field * fieldDat,double ** field, int i, int j) {
//     double L[3];
//     double x[4];
//     int inc = 1;
//
//     if (i > fieldDat->fieldLatLen / 2) inc = -inc;
//
//     x[0] = fieldDat->lat[i];
//     x[1] = fieldDat->lat[i + inc];
//     x[2] = fieldDat->lat[i + inc*2];
//     x[3] = fieldDat->lat[i + inc*3];
//
//     //solve for each quadratic term
//     L[0] = (x[0] - x[2])*(x[0] - x[3]) / ((x[1] - x[2])*(x[1] - x[3])) * field[i+inc][j];             //L1
//     L[1] = (x[0] - x[1])*(x[0] - x[3]) / ((x[2] - x[1])*(x[2] - x[3])) * field[i + inc*2][j];             //L2
//     L[2] = (x[0] - x[1])*(x[0] - x[2]) / ((x[3] - x[1])*(x[3] - x[2]))* field[i + inc * 3][j];             //L3
//
//     return L[0]+L[1]+L[2];
// };
//
// double linearInterp4(Field * field, int i, int j);
// inline double linearInterp4Array(Field * fieldDat, double ** field, int i, int j) __attribute__((always_inline));
// inline double linearInterp4Array(Field * fieldDat, double ** field, int i, int j) {
//     //fourth order accurate forward difference method
//     if (i < fieldDat->fieldLatLen / 2) {
//         return field[i + 1][j] + (25. / 12.* field[i + 1][j] - 4 * field[i + 2][j] + 3 * field[i + 3][j] - 4 / 3. * field[i + 4][j] + .25 * field[i + 5][j]);
//     }
//     else {
//         return field[i - 1][j] + (25. / 12.* field[i - 1][j] - 4 * field[i - 2][j] + 3 * field[i - 3][j] - 4 / 3. * field[i - 4][j] + .25 * field[i - 5][j]);
//     }
// };
//
//
//
// double lagrangeInterp2(Field * field, int i, int j);
// double lagrangeInterp3(Field * field, int i, int j);
// // double lagrangeInterp3Array(Field * fieldDat, double **, int i, int j);
// //double multilinFit(Field * field, int i, int j);
// inline double lagrangeInterp3Array(Field * fieldDat, double ** field, int i, int j) __attribute__((always_inline));
// inline double lagrangeInterp3Array(Field * fieldDat, double ** field, int i, int j)  {
//     double L[4];
//     double x[5];
//     int inc = 1;
//
//     if (i > fieldDat->fieldLatLen / 2) inc = -inc;
//
//
//     x[0] = fieldDat->lat[i];
//     x[1] = fieldDat->lat[i + inc];
//     x[2] = fieldDat->lat[i + inc*2];
//     x[3] = fieldDat->lat[i + inc*3];
//     x[4] = fieldDat->lat[i + inc*4];
//
//     //solve for each quadratic term
//     L[0] = (x[0] - x[2])*(x[0] - x[3])*(x[0] - x[4]) / ((x[1] - x[2])*(x[1] - x[3])*(x[1] - x[4])) * field[i + inc][j];             //L1
//     L[1] = (x[0] - x[1])*(x[0] - x[3])*(x[0] - x[4]) / ((x[2] - x[1])*(x[2] - x[3])*(x[2] - x[4])) * field[i + inc*2][j];             //L2
//     L[2] = (x[0] - x[1])*(x[0] - x[2])*(x[0] - x[4]) / ((x[3] - x[1])*(x[3] - x[2])*(x[3] - x[4])) * field[i + inc*3][j];             //L3
//     L[3] = (x[0] - x[1])*(x[0] - x[2])*(x[0] - x[3]) / ((x[4] - x[1])*(x[4] - x[2])*(x[4] - x[3])) * field[i + inc*4][j];
//
//     return L[0] + L[1] + L[2] + L[3];
// };
//
// inline double lagrangeInterp3ArrayCenter(Field * fieldDat, double ** field, int i, int j) __attribute__((always_inline));
// inline double lagrangeInterp3ArrayCenter(Field * fieldDat, double ** field, int i, int j)  {
//     double L[4];
//     double x[5];
//     int inc = 1;
//
//     if (i > fieldDat->fieldLatLen / 2) inc = -inc;
//
//
//     x[0] = fieldDat->lat[i];
//     x[1] = x[0]-fieldDat->dLat*2*inc;
//     x[2] = x[0]-fieldDat->dLat*inc;
//     x[3] = x[0]+fieldDat->dLat*inc;
//     x[4] = x[0]+fieldDat->dLat*2*inc;
//
//     //solve for each quadratic term
//     L[0] = (x[0] - x[2])*(x[0] - x[3])*(x[0] - x[4]) / ((x[1] - x[2])*(x[1] - x[3])*(x[1] - x[4])) * field[i+inc*2][j];             //L1
//     L[1] = (x[0] - x[1])*(x[0] - x[3])*(x[0] - x[4]) / ((x[2] - x[1])*(x[2] - x[3])*(x[2] - x[4])) * field[i+inc][j];             //L2
//
//     if (j<(fieldDat->fieldLonLen/2)) {
//         L[2] = (x[0] - x[1])*(x[0] - x[2])*(x[0] - x[4]) / ((x[3] - x[1])*(x[3] - x[2])*(x[3] - x[4])) * field[i+inc][j+(fieldDat->fieldLonLen / 2)];                         //L3
//         L[3] = (x[0] - x[1])*(x[0] - x[2])*(x[0] - x[3]) / ((x[4] - x[1])*(x[4] - x[2])*(x[4] - x[3])) * field[i+inc*2][j+(fieldDat->fieldLonLen / 2)];
//     }
//     else {
//         L[2] = (x[0] - x[1])*(x[0] - x[2])*(x[0] - x[4]) / ((x[3] - x[1])*(x[3] - x[2])*(x[3] - x[4])) * field[i+inc][j-(fieldDat->fieldLonLen / 2)];                         //L3
//         L[3] = (x[0] - x[1])*(x[0] - x[2])*(x[0] - x[3]) / ((x[4] - x[1])*(x[4] - x[2])*(x[4] - x[3])) * field[i+inc*2][j-(fieldDat->fieldLonLen / 2)];
//     }
//
//     return L[0] + L[1] + L[2] + L[3];
// };

#endif
