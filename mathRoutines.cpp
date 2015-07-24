#include "mathRoutines.h"
#include "field.h"
#include "globals.h"
#include "solver.h"
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_multifit.h>

double linearInterp1(Field * field, int i, int j) {
	
	//third order accurate forward difference method
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

	//third order accurate forward difference method
	if (i < field->fieldLatLen / 2) {
		return field->solution[i + 1][j] + (1.5* field->solution[i + 1][j] - 2 * field->solution[i + 2][j] + 0.5* field->solution[i + 3][j]);
	}
	else {
		return field->solution[i - 1][j] + (1.5* field->solution[i - 1][j] - 2 * field->solution[i - 2][j] + 0.5* field->solution[i - 3][j]);

	}

	

	//return 0;
};

double linearInterp4(Field * field, int i, int j) {
	//third order accurate forward difference method
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


//double multilinFit(Field * field, int i, int j) {
//	
//	
//	int lon1;
//	int lon2;
//	int lon3;
//
//	double v1;
//	double v2;
//	double v3;
//	double v4;
//	double v5;
//	double v6;
//
//
//	int lat1, lat2;
//
//	if (i < field->fieldLatLen / 2) {
//		lat1 = field->lat[i + 1];
//		lat2 = field->lat[i + 2];
//
//		if (j == 0) {
//			lon1 = field->lon[j] - field->dLon; // field->lon[field->fieldLonLen - 1];
//			lon2 = field->lon[j];
//			lon3 = field->lon[j + 1];
//
//			v1 = field->solution[i + 1][field->fieldLonLen - 1];
//			v2 = field->solution[i + 1][j];
//			v3 = field->solution[i + 1][j + 1];
//			v4 = field->solution[i + 2][field->fieldLonLen - 1];
//			v5 = field->solution[i + 2][j];
//			v6 = field->solution[i + 2][j + 1];
//		}
//
//		else if (j == field->fieldLonLen - 1) {
//			lon1 = field->lon[j - 1];
//			lon2 = field->lon[j];
//			lon3 = field->lon[j] + field->dLon;//field->lon[0];
//
//			v1 = field->solution[i + 1][j - 1];
//			v2 = field->solution[i + 1][j];
//			v3 = field->solution[i + 1][0];
//			v4 = field->solution[i + 2][j - 1];
//			v5 = field->solution[i + 2][j];
//			v6 = field->solution[i + 2][0];
//		}
//		else {
//			lon1 = field->lon[j - 1];
//			lon2 = field->lon[j];
//			lon3 = field->lon[j + 1];
//
//			v1 = field->solution[i + 1][j - 1];
//			v2 = field->solution[i + 1][j];
//			v3 = field->solution[i + 1][j + 1];
//			v4 = field->solution[i + 2][j - 1];
//			v5 = field->solution[i + 2][j];
//			v6 = field->solution[i + 2][j + 1];
//		}
//	}
//	else {
//		lat1 = field->lat[i - 2];
//		lat2 = field->lat[i - 1];
//
//		if (j == 0) {
//			lon1 = field->lon[j] - field->dLon; // field->lon[field->fieldLonLen - 1];
//			lon2 = field->lon[j];
//			lon3 = field->lon[j + 1];
//
//			v1 = field->solution[i - 2][field->fieldLonLen - 1];
//			v2 = field->solution[i - 2][j];
//			v3 = field->solution[i - 2][j + 1];
//			v4 = field->solution[i - 1][field->fieldLonLen - 1];
//			v5 = field->solution[i - 1][j];
//			v6 = field->solution[i - 1][j + 1];
//		}
//
//		else if (j == field->fieldLonLen - 1) {
//			lon1 = field->lon[j - 1];
//			lon2 = field->lon[j];
//			lon3 = field->lon[j] + field->dLon;//field->lon[0];
//
//			v1 = field->solution[i - 2][j - 1];
//			v2 = field->solution[i - 2][j];
//			v3 = field->solution[i - 2][0];
//			v4 = field->solution[i - 1][j - 1];
//			v5 = field->solution[i - 1][j];
//			v6 = field->solution[i - 1][0];
//		}
//		else {
//			lon1 = field->lon[j - 1];
//			lon2 = field->lon[j];
//			lon3 = field->lon[j + 1];
//
//			v1 = field->solution[i - 2][j - 1];
//			v2 = field->solution[i - 2][j];
//			v3 = field->solution[i - 2][j + 1];
//			v4 = field->solution[i - 1][j - 1];
//			v5 = field->solution[i - 1][j];
//			v6 = field->solution[i - 1][j + 1];
//		}
//	}
//
//
//	gsl_matrix * X = gsl_matrix_alloc(4, 3);
//
//	gsl_matrix_set(X, 0, 0, 1);
//	gsl_matrix_set(X, 1, 0, 1);
//	gsl_matrix_set(X, 2, 0, 1);
//	gsl_matrix_set(X, 3, 0, 1);
//	//gsl_matrix_set(X, 4, 0, 1);
//	//gsl_matrix_set(X, 5, 0, 1);
//
//	//column 2
//	gsl_matrix_set(X, 0, 1, lat1);
//	gsl_matrix_set(X, 1, 1, lat1);
//	//gsl_matrix_set(X, 2, 1, lat1);
//	gsl_matrix_set(X, 2, 1, lat2);
//	gsl_matrix_set(X, 3, 1, lat2);
//	//gsl_matrix_set(X, 5, 1, lat2);
//
//	//column 2
//	gsl_matrix_set(X, 0, 2, lon1);
//	gsl_matrix_set(X, 1, 2, lon3);
//	gsl_matrix_set(X, 2, 2, lon1);
//	gsl_matrix_set(X, 3, 2, lon3);
//	//gsl_matrix_set(X, 4, 2, lon2);
//	//gsl_matrix_set(X, 5, 2, lon3);
//
//	//Matrix shape:
//	// 1 lat1 lon1
//	// 1 lat1 lon2
//	// 1 lat1 lon3
//	// 1 lat2 lon1
//	// 1 lat2 lon2
//	// 1 lat2 lon3
//
//	gsl_vector * v = gsl_vector_alloc(4);
//
//	gsl_vector_set(v, 0, v1);
//	//gsl_vector_set(v, 1, v2);
//	gsl_vector_set(v, 1, v3);
//	gsl_vector_set(v, 2, v4);
//	//gsl_vector_set(v, 4, v5);
//	gsl_vector_set(v, 3, v6);
//
//	gsl_vector * c = gsl_vector_alloc(3);
//
//	gsl_vector_set(c, 0, 0);
//	gsl_vector_set(c, 1, 0);
//	gsl_vector_set(c, 2, 0);
//
//	//Create workspace
//
//	gsl_matrix * cov = gsl_matrix_alloc(3, 3);
//
//	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc(4, 3);
//
//	double chisq;
//
//	int sumSq;
//
//	sumSq = gsl_multifit_linear(X, v, c, cov, &chisq, work);
//
//	double vnew = gsl_vector_get(c, 0) + gsl_vector_get(c, 1) * field->lat[i] + gsl_vector_get(c, 2) * field->lon[j];
//
//	gsl_multifit_linear_free(work);
//	gsl_vector_free(c);
//	gsl_vector_free(v);
//	gsl_matrix_free(X);
//	gsl_matrix_free(cov);
//
//	return vnew;
//}
