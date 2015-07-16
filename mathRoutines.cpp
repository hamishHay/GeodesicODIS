#include "mathRoutines.h"
#include "field.h"
#include "globals.h"

using namespace std;

double linearInterp(Field * field, int i, int j) {
	//third order accurate forward difference method
	if (i < field->fieldLatLen / 2) {
		return -(11. / 6.* field->solution[i + 1][j] - 3 * field->solution[i + 2][j] + 1.5* field->solution[i + 3][j] - 1 / 3. * field->solution[i + 4][j]);
	}
	else {
		return -(11. / 6.* field->solution[i - 1][j] - 3 * field->solution[i - 2][j] + 1.5* field->solution[i - 3][j] - 1 / 3. * field->solution[i - 4][j]);
	}
}

double lagrangeInterp(Field * field, int i, int j) {
	double L[3];
	double P2 = 0.;
	double x[4];
	int inc = 3;
	int count = 0;
	
	if (i < field->fieldLatLen / 2) {
		inc = 1*inc;
	}
	else {
		inc = -inc;
	}
	//define x positions (i.e. colat) in radians
	count = 1;
	x[0] = (90 - field->lat[i]);// *radConv;
	for (int n = inc/abs(inc); abs(n) <= abs(inc*3); n+=inc) {
		x[count] = (90 - field->lat[i + n]);// *radConv;
		count++;
	}

	//solve for each quadratic term
	L[0] = (x[0] - x[2])*(x[0] - x[3]) / ((x[1] - x[2])*(x[1] - x[3])); //L1
	L[1] = (x[0] - x[1])*(x[0] - x[3]) / ((x[2] - x[1])*(x[2] - x[3])); //L2
	L[2] = (x[0] - x[1])*(x[0] - x[2]) / ((x[3] - x[1])*(x[3] - x[2])); //L3

	//interpolated value calculated via summation of each L term.
	count = 0;
	for (int n = inc/abs(inc); abs(n) <= abs(inc * 3); n += inc){
		P2 += field->solution[i + n][j] * L[count];
		count++;
	}

	return P2;
}
