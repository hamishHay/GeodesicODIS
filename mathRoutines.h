#ifndef MATHROUTINES_H
#define MATHROUTINES_H

#include "field.h"
#include "solver.h"
#include <math.h>

double arcLength(double lat1, double lat2, double lon1, double lon2);

inline double gamm1n(double xx);
inline double gamm1n(double xx) {
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,
                          -86.50532032941677,
                          24.01409824083091,
                          -1.231739572450155,
                          0.1208650973866179e-2,
                          -0.5395239384953e-5};
    int j;

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0; j<=5; j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

inline double factrl(int n)
{
    double gammln(float xx);

    static int ntop=4;
    static float a[33]={1.0,
                        1.0,
                        2.0,
                        6.0,
                        24.0,
                        120.0,
                        720.0,
                        5040.0,
                        40320.0,
                        362880.0,
                        3628800.0,
                        39916800.0,
                        479001600.0,
                        6227020800.0,
                        87178291200.0,
                        1307674368000.0,
                        20922789888000.0,
                        355687428096000.0,
                        6402373705728000.0,
                        1.21645100408832e+17,
                        2.43290200817664e+18,
                        5.109094217170944e+19,
                        1.1240007277776077e+21,
                        2.585201673888498e+22,
                        6.204484017332394e+23,
                        1.5511210043330984e+25,
                        4.032914611266057e+26,
                        1.0888869450418352e+28,
                        3.048883446117138e+29,
                        8.841761993739701e+30,
                        2.652528598121911e+32,
                        8.222838654177924e+33,
                        2.6313083693369355e+35};

    int j;

    if (n > 32) return exp(gamm1n(n+1.0));
    // while (ntop<n) {
    //   j=ntop++;
    //   a[ntop]=a[j]*ntop;
    // }
    return a[n];
}

inline double legendreNorm(int l, int m);
inline double legendreNorm(int l, int m) {
    double norm[9][9] = {{1.0, 0, 0, 0, 0, 0, 0, 0, 0},
                         {1.73205080757, 1.73205080757, 0, 0, 0, 0, 0, 0, 0},
                         {2.2360679775, 1.29099444874, 0.645497224368, 0, 0, 0, 0, 0, 0},
                         {2.64575131106, 1.08012344973, 0.341565025532, 0.139443337756, 0, 0, 0, 0, 0},
                         {3.0, 0.948683298051, 0.22360679775, 0.0597614304667, 0.0211288563682, 0, 0, 0, 0},
                         {3.31662479036, 0.856348838578, 0.161834718743, 0.0330343736322, 0.00778627653585, 0.00246223683452, 0, 0, 0},
                         {3.60555127546, 0.786795792469, 0.124403337882, 0.020733889647, 0.0037854730215, 0.000807065559928, 0.000232979759139, 0, 0},
                         {3.87298334621, 0.731925054711, 0.0996023841112, 0.0140859042455, 0.00212352996432, 0.000353921660721, 6.94097482422e-05, 1.8550535516e-05, 0},
                         {4.12310562562, 0.687184270936, 0.0821342300508, 0.0101100248374, 0.00130519859416, 0.000180998479074, 2.79286716592e-05, 5.09905448963e-06, 1.27476362241e-06}};

    return norm[l][m];
}

// Numerical recipes function for calculating the associated Legendre
// function of x.
inline double assLegendre(int l, int m, double x) __attribute__((always_inline));
inline double assLegendre(int l, int m, double x) {
    double fact, pll, pmm, pmmp1, somx2, normalise;
    double delta = 0.0;
    int i, ll;

    if (m==0) delta = 1.0;

    // normalise = legendreNorm(l,m);
    normalise =  sqrt((2.0-delta)*(2.0*(double)l+1.0)*factrl(l-m)/factrl(l+m));

    // if (m%2 != 0) normalise = -normalise;

    // std::cout<<"l="<<l<<" m="<<m<<" norm="<<normalise<<std::endl;

    pmm = 1.0;
    if (m > 0)
    {
        somx2 = sqrt((1.0 - x)*(1.0+x));
        fact = 1.0;
        for (i=1; i<=m; i++)
        {
            pmm *= -fact*somx2;
            fact += 2.0;
        }
    }
    if (l == m) {
        return normalise*pmm;
    }
    else {
        pmmp1 = x*(2*m+1)*pmm;
        if (l == (m+1)) {
            return normalise*pmmp1;
        }
        else {
            for (ll=m+2; ll<=l; ll++) {
                pll = (x*(2*ll-1)*pmmp1 - (ll+m-1)*pmm)/(ll-m);
                pmm = pmmp1;
                pmmp1 = pll;
            }
            return normalise*pll;
        }
    }
};

inline double fastSin(double x);
inline double fastSin(double x) {
    double ans;
    double xxx = x*x*x;
    double xxxxx = xxx*x*x;

    ans = x;
    ans -= (xxx)/factrl(3);
    ans += (xxxxx)/factrl(5);
}

//func returns value n at given i and j based on Lagrange Interpolation
//of points n+1, n+2, and n+3 for a 2nd degree polynomial
inline double linearInterp1Array(Field * field, double ** solution, int i, int j) __attribute__((always_inline));
inline double linearInterp1Array(Field * field, double ** solution, int i, int j) {
    //first order accurate forward difference method
    if (i < field->fieldLatLen / 2) {
        return solution[i + 1][j] + (solution[i + 1][j] - solution[i + 2][j]);
    }
    else {
        return solution[i - 1][j] + (solution[i - 1][j] - solution[i - 2][j]);
    }
};





double linearInterp2(Field * field, int i, int j);

inline double lagrangeInterp2Array(Field * fieldDat, double ** field, int i, int j) __attribute__((always_inline));
inline double lagrangeInterp2Array(Field * fieldDat,double ** field, int i, int j) {
    double L[3];
    double x[4];
    int inc = 1;

    if (i > fieldDat->fieldLatLen / 2) inc = -inc;

    x[0] = fieldDat->lat[i];
    x[1] = fieldDat->lat[i + inc];
    x[2] = fieldDat->lat[i + inc*2];
    x[3] = fieldDat->lat[i + inc*3];

    //solve for each quadratic term
    L[0] = (x[0] - x[2])*(x[0] - x[3]) / ((x[1] - x[2])*(x[1] - x[3])) * field[i+inc][j];             //L1
    L[1] = (x[0] - x[1])*(x[0] - x[3]) / ((x[2] - x[1])*(x[2] - x[3])) * field[i + inc*2][j];             //L2
    L[2] = (x[0] - x[1])*(x[0] - x[2]) / ((x[3] - x[1])*(x[3] - x[2]))* field[i + inc * 3][j];             //L3

    return L[0]+L[1]+L[2];
};

double linearInterp4(Field * field, int i, int j);
inline double linearInterp4Array(Field * fieldDat, double ** field, int i, int j) __attribute__((always_inline));
inline double linearInterp4Array(Field * fieldDat, double ** field, int i, int j) {
    //fourth order accurate forward difference method
    if (i < fieldDat->fieldLatLen / 2) {
        return field[i + 1][j] + (25. / 12.* field[i + 1][j] - 4 * field[i + 2][j] + 3 * field[i + 3][j] - 4 / 3. * field[i + 4][j] + .25 * field[i + 5][j]);
    }
    else {
        return field[i - 1][j] + (25. / 12.* field[i - 1][j] - 4 * field[i - 2][j] + 3 * field[i - 3][j] - 4 / 3. * field[i - 4][j] + .25 * field[i - 5][j]);
    }
};



double lagrangeInterp2(Field * field, int i, int j);
double lagrangeInterp3(Field * field, int i, int j);
// double lagrangeInterp3Array(Field * fieldDat, double **, int i, int j);
//double multilinFit(Field * field, int i, int j);
inline double lagrangeInterp3Array(Field * fieldDat, double ** field, int i, int j) __attribute__((always_inline));
inline double lagrangeInterp3Array(Field * fieldDat, double ** field, int i, int j)  {
    double L[4];
    double x[5];
    int inc = 1;

    if (i > fieldDat->fieldLatLen / 2) inc = -inc;


    x[0] = fieldDat->lat[i];
    x[1] = fieldDat->lat[i + inc];
    x[2] = fieldDat->lat[i + inc*2];
    x[3] = fieldDat->lat[i + inc*3];
    x[4] = fieldDat->lat[i + inc*4];

    //solve for each quadratic term
    L[0] = (x[0] - x[2])*(x[0] - x[3])*(x[0] - x[4]) / ((x[1] - x[2])*(x[1] - x[3])*(x[1] - x[4])) * field[i + inc][j];             //L1
    L[1] = (x[0] - x[1])*(x[0] - x[3])*(x[0] - x[4]) / ((x[2] - x[1])*(x[2] - x[3])*(x[2] - x[4])) * field[i + inc*2][j];             //L2
    L[2] = (x[0] - x[1])*(x[0] - x[2])*(x[0] - x[4]) / ((x[3] - x[1])*(x[3] - x[2])*(x[3] - x[4])) * field[i + inc*3][j];             //L3
    L[3] = (x[0] - x[1])*(x[0] - x[2])*(x[0] - x[3]) / ((x[4] - x[1])*(x[4] - x[2])*(x[4] - x[3])) * field[i + inc*4][j];

    return L[0] + L[1] + L[2] + L[3];
};

inline double lagrangeInterp3ArrayCenter(Field * fieldDat, double ** field, int i, int j) __attribute__((always_inline));
inline double lagrangeInterp3ArrayCenter(Field * fieldDat, double ** field, int i, int j)  {
    double L[4];
    double x[5];
    int inc = 1;

    if (i > fieldDat->fieldLatLen / 2) inc = -inc;


    x[0] = fieldDat->lat[i];
    x[1] = x[0]-fieldDat->dLat*2*inc;
    x[2] = x[0]-fieldDat->dLat*inc;
    x[3] = x[0]+fieldDat->dLat*inc;
    x[4] = x[0]+fieldDat->dLat*2*inc;

    //solve for each quadratic term
    L[0] = (x[0] - x[2])*(x[0] - x[3])*(x[0] - x[4]) / ((x[1] - x[2])*(x[1] - x[3])*(x[1] - x[4])) * field[i+inc*2][j];             //L1
    L[1] = (x[0] - x[1])*(x[0] - x[3])*(x[0] - x[4]) / ((x[2] - x[1])*(x[2] - x[3])*(x[2] - x[4])) * field[i+inc][j];             //L2

    if (j<(fieldDat->fieldLonLen/2)) {
        L[2] = (x[0] - x[1])*(x[0] - x[2])*(x[0] - x[4]) / ((x[3] - x[1])*(x[3] - x[2])*(x[3] - x[4])) * field[i+inc][j+(fieldDat->fieldLonLen / 2)];                         //L3
        L[3] = (x[0] - x[1])*(x[0] - x[2])*(x[0] - x[3]) / ((x[4] - x[1])*(x[4] - x[2])*(x[4] - x[3])) * field[i+inc*2][j+(fieldDat->fieldLonLen / 2)];
    }
    else {
        L[2] = (x[0] - x[1])*(x[0] - x[2])*(x[0] - x[4]) / ((x[3] - x[1])*(x[3] - x[2])*(x[3] - x[4])) * field[i+inc][j-(fieldDat->fieldLonLen / 2)];                         //L3
        L[3] = (x[0] - x[1])*(x[0] - x[2])*(x[0] - x[3]) / ((x[4] - x[1])*(x[4] - x[2])*(x[4] - x[3])) * field[i+inc*2][j-(fieldDat->fieldLonLen / 2)];
    }

    return L[0] + L[1] + L[2] + L[3];
};

#endif
