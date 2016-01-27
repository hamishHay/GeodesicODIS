#ifndef MATHROUTINES_H
#define MATHROUTINES_H

#include "field.h"
#include "solver.h"

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
