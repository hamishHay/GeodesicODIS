#ifndef MATHROUTINES_H
#define MATHROUTINES_H

#include "field.h"
#include "solver.h"

//func returns value n at given i and j based on Lagrange Interpolation
//of points n+1, n+2, and n+3 for a 2nd degree polynomial
double linearInterp1(Field * field, int i, int j);
double linearInterp2(Field * field, int i, int j);
double linearInterp4(Field * field, int i, int j);
double lagrangeInterp(Field * field, int i, int j);
double lagrangeInterp2(Field * field, int i, int j);
//double multilinFit(Field * field, int i, int j);

#endif
