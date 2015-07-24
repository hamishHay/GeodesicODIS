#include "mathRoutines.h"
#include "field.h"
#include "globals.h"
#include "solver.h"

double linearInterp1(Field * field, int i, int j) {
	
	//first order accurate forward difference method
	if (i < field->fieldLatLen / 2) {
		return field->solution[i + 1][j] + (field->solution[i + 1][j] - field->solution[i + 2][j]);
		//return field->solution[i + 1][j] + (1.5* field->solution[i + 1][j] - 2 * field->solution[i + 2][j] + 0.5* field->solution[i + 3][j]);
		//return field->solution[i + 1][j] + (11. / 6.* field->solution[i + 1][j] - 3 * field->solution[i + 2][j] + 1.5* field->solution[i + 3][j] - 1 / 3. * field->solution[i + 4][j]);
	}
	else {
		return field->solution[i - 1][j] + (field->solution[i - 1][j] - field->solution[i - 2][j]);
		//return field->solution[i - 1][j] + (1.5* field->solution[i - 1][j] - 2 * field->solution[i - 2][j] + 0.5* field->solution[i - 3][j]);
		//return field->solution[i - 1][j] + (11. / 6.* field->solution[i - 1][j] - 3 * field->solution[i - 2][j] + 1.5* field->solution[i - 3][j] - 1 / 3. * field->solution[i - 4][j]);
	}
};

double linearInterp2(Field * field, int i, int j) {

	//second order accurate forward difference method
	if (i < field->fieldLatLen / 2) {
		return field->solution[i + 1][j] + (1.5* field->solution[i + 1][j] - 2 * field->solution[i + 2][j] + 0.5* field->solution[i + 3][j]);
	}
	else {
		return field->solution[i - 1][j] + (1.5* field->solution[i - 1][j] - 2 * field->solution[i - 2][j] + 0.5* field->solution[i - 3][j]);

	}
};

double linearInterp4(Field * field, int i, int j) {
	//fourth order accurate forward difference method
	if (i < field->fieldLatLen / 2) {
		return field->solution[i + 1][j] + (25. / 12.* field->solution[i + 1][j] - 4 * field->solution[i + 2][j] + 3 * field->solution[i + 3][j] - 4 / 3. * field->solution[i + 4][j] + .25 * field->solution[i + 5][j]);
	}
	else {
		return field->solution[i - 1][j] + (25. / 12.* field->solution[i - 1][j] - 4 * field->solution[i - 2][j] + 3 * field->solution[i - 3][j] - 4 / 3. * field->solution[i - 4][j] + .25 * field->solution[i - 5][j]);
		//return field->solution[i - 1][j] + (1.5* field->solution[i - 1][j] - 2 * field->solution[i - 2][j] + 0.5* field->solution[i - 3][j]);
		//return field->solution[i - 1][j] + (11. / 6.* field->solution[i - 1][j] - 3 * field->solution[i - 2][j] + 1.5* field->solution[i - 3][j] - 1 / 3. * field->solution[i - 4][j]);
	}
};

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
};

double lagrangeInterp2(Field * field, int i, int j) {
	double L[3];
	double x[4];
	int inc = 1;

	if (i > field->fieldLatLen / 2) inc = -inc;

	x[0] = field->lat[i];
	x[1] = field->lat[i + inc];
	x[2] = field->lat[i + inc*2];
	x[3] = field->lat[i + inc*3];

	//solve for each quadratic term
	L[0] = (x[0] - x[2])*(x[0] - x[3]) / ((x[1] - x[2])*(x[1] - x[3])) * field->solution[i+inc][j]; //L1
	L[1] = (x[0] - x[1])*(x[0] - x[3]) / ((x[2] - x[1])*(x[2] - x[3])) * field->solution[i + inc*2][j]; //L2
	L[2] = (x[0] - x[1])*(x[0] - x[2]) / ((x[3] - x[1])*(x[3] - x[2]))* field->solution[i + inc * 3][j]; //L3

	return L[0]+L[1]+L[2];
};