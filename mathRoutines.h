#ifndef MATHROUTINES_H
#define MATHROUTINES_H

#include "field.h"
#include "solver.h"
#include <math.h>

inline int factorial(int n) __attribute__((always_inline));
inline int factorial(int n)
{
    if(n > 1)
        return n * factorial(n - 1);
    else
        return 1;
}

// Numerical recipes function for calculating the associated Legendre
// function of x.
inline double assLegendre(int l, int m, double x) __attribute__((always_inline));
inline double assLegendre(int l, int m, double x) {
  double fact, pll, pmm, pmmp1, somx2, normalise;
  double delta = 0.0;
  int i, ll;

  if (m==0) delta = 1.0;

  normalise = sqrt((2.0-delta)*(2.0*l+1.0)*factorial(l-m)/factorial(l+m));

  pmm = 1.0;
  if (m > 0)
  {
    somx2 = sqrt((1.0 - x)*(1.0+x));
    fact = 1.0;
    for (i=1; i<m; i++)
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
      for (ll=m+2; ll<=1; ll++) {
        pll = (x*(2*ll-1)*pmmp1 - (ll+m-1)*pmm)/(ll-m);
        pmm = pmmp1;
        pmmp1 = pll;
      }
      return normalise*pll;
    }
  }
};

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
	L[0] = (x[0] - x[2])*(x[0] - x[3]) / ((x[1] - x[2])*(x[1] - x[3])) * field[i+inc][j]; //L1
	L[1] = (x[0] - x[1])*(x[0] - x[3]) / ((x[2] - x[1])*(x[2] - x[3])) * field[i + inc*2][j]; //L2
	L[2] = (x[0] - x[1])*(x[0] - x[2]) / ((x[3] - x[1])*(x[3] - x[2]))* field[i + inc * 3][j]; //L3

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
	L[0] = (x[0] - x[2])*(x[0] - x[3])*(x[0] - x[4]) / ((x[1] - x[2])*(x[1] - x[3])*(x[1] - x[4])) * field[i + inc][j]; //L1
	L[1] = (x[0] - x[1])*(x[0] - x[3])*(x[0] - x[4]) / ((x[2] - x[1])*(x[2] - x[3])*(x[2] - x[4])) * field[i + inc*2][j]; //L2
	L[2] = (x[0] - x[1])*(x[0] - x[2])*(x[0] - x[4]) / ((x[3] - x[1])*(x[3] - x[2])*(x[3] - x[4])) * field[i + inc*3][j]; //L3
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
	L[0] = (x[0] - x[2])*(x[0] - x[3])*(x[0] - x[4]) / ((x[1] - x[2])*(x[1] - x[3])*(x[1] - x[4])) * field[i+inc*2][j]; //L1
	L[1] = (x[0] - x[1])*(x[0] - x[3])*(x[0] - x[4]) / ((x[2] - x[1])*(x[2] - x[3])*(x[2] - x[4])) * field[i+inc][j]; //L2

  if (j<(fieldDat->fieldLonLen/2)) {
    L[2] = (x[0] - x[1])*(x[0] - x[2])*(x[0] - x[4]) / ((x[3] - x[1])*(x[3] - x[2])*(x[3] - x[4])) * field[i+inc][j+(fieldDat->fieldLonLen / 2)]; //L3
  	L[3] = (x[0] - x[1])*(x[0] - x[2])*(x[0] - x[3]) / ((x[4] - x[1])*(x[4] - x[2])*(x[4] - x[3])) * field[i+inc*2][j+(fieldDat->fieldLonLen / 2)];
  }
  else {
    L[2] = (x[0] - x[1])*(x[0] - x[2])*(x[0] - x[4]) / ((x[3] - x[1])*(x[3] - x[2])*(x[3] - x[4])) * field[i+inc][j-(fieldDat->fieldLonLen / 2)]; //L3
  	L[3] = (x[0] - x[1])*(x[0] - x[2])*(x[0] - x[3]) / ((x[4] - x[1])*(x[4] - x[2])*(x[4] - x[3])) * field[i+inc*2][j-(fieldDat->fieldLonLen / 2)];
  }

	return L[0] + L[1] + L[2] + L[3];
};

#endif
